#include "stp_model_solver.h"
//#include "stp_benders.h"
#include "stp_ratio_heuristic.h"



void assignment_based_model(StpInstance &xIns) {
	clock_t startTime, endTime;
	startTime = clock();
	//Define variable x_ik: (=1) if job i is assigned to batch k
	GRBVar** x = new GRBVar * [STP_N + 1];
	for (int i = 1; i <= STP_N; ++i) {
		x[i] = new GRBVar[STP_N + 1];
	}
	// Define variable y_i: the probability of component i to be executed 
	GRBVar* y = new GRBVar[STP_N + 1];

	//Define variable z_ik
	GRBVar** z = nullptr;
	z = new GRBVar * [STP_N + 1];
	for (int i = 0; i <= STP_N; ++i) {
		z[i] = new GRBVar[STP_N + 1];
	}

	// Define variable u_k: (=1) if the k-th batch is used
	GRBVar* u = new GRBVar[STP_N + 1];

	GRBVar* f = new GRBVar[STP_N + 1];
	try {
		GRBEnv env = GRBEnv();
		GRBModel model = GRBModel(env);
		
		//Define variable x_ik: (=1) if job i is assigned to batch k
		for (int i = 1; i <= STP_N; ++i) {
			for (int k = 1; k <= STP_N; ++k) {
				ostringstream vname;
				vname << "x_" << i << "_" << k;
				x[i][k] = model.addVar(0, 1, 0, GRB_BINARY,vname.str());
			}
		}
		// Define variable y_i: the probability of component i to be executed 
		for (int i = 1; i <= STP_N; ++i) {
			ostringstream vname;
			vname << "y_" << i;
			y[i] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, vname.str());
		}
		//Define variable z_ik
		for (int i = 0; i <= STP_N; ++i) {
			for (int k = 1; k <= STP_N; ++k) {
				ostringstream vname;
				vname << "z_" << i << "_" << k;
				z[i][k] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, vname.str());
			}
		}
		// Define variable u_k: (=1) if the k-th batch is used
		for (int k = 1; k <= STP_N; ++k) {
			ostringstream vname;
			vname << "u_" << k;
			u[k] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, vname.str());
		}

		for (int k = 1; k <= STP_N; ++k) {
			ostringstream vname;
			vname << "f_" << k;
			f[k] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, vname.str());
		}
		//////////////////////////////////////////
		//Add constraints
		//Component assignment
		for (int i = 1; i <= STP_N; i++){
			GRBLinExpr expr = 0;
			for (int k = 1; k <= STP_N; ++k) {
				expr += x[i][k];
			}
			model.addConstr(expr, GRB_EQUAL, 1);
		}

		//Resource constraints
		for (int k = 1; k <= STP_N; k++){
			GRBLinExpr expr = 0;
			for (int i = 1; i <= STP_N; i++){
				expr += x[i][k] * xIns.mComp[i].mRes;
			}
			model.addConstr(expr, GRB_LESS_EQUAL, STP_RESCAP);
		}

		for (int i = 1; i <= STP_N; i++){
			for (int k = 1; k <= STP_N; k++){
				//GRBLinExpr expr = 0;
				//for (int v = 1; v <=k; v++){
					//expr += x[i][v];
				//}
				model.addConstr(y[i], GRB_GREATER_EQUAL, z[STP_N][k] - 1 + x[i][k]);
			}
		}

		for (int i = 0; i <= STP_N; i++){
			model.addConstr(z[i][1], GRB_EQUAL, 1);
		}
		for (int k = 2; k <= STP_N; k++){
			model.addConstr(z[0][k], GRB_EQUAL, z[STP_N][k - 1]);
		}

		for (int i = 1; i <= STP_N; i++){
			for (int k = 2; k <= STP_N; k++){
			model.addConstr(z[i][k], GRB_GREATER_EQUAL, z[i - 1][k]-x[i][k-1]);
			model.addConstr(z[i][k], GRB_GREATER_EQUAL, xIns.mComp[i].mProb* z[i - 1][k]);
			}
		}

		for (int k = 1; k <= STP_N; k++){
			for (int i = 1; i <= STP_N; i++){
				model.addConstr(u[k], GRB_GREATER_EQUAL, z[STP_N][k] - 1 + x[i][k]);
			}
		}
		
		for (int k = 1; k <= STP_N; k++) {
			GRBLinExpr exp = 0;
			for (int i = 1; i <= STP_N; i++) {
				exp += x[i][k];
				model.addConstr(x[i][k], GRB_LESS_EQUAL, f[k]);
			}
			model.addConstr(exp, GRB_GREATER_EQUAL, f[k]);
		}
		for (int k = 2; k <= STP_N; k++) {
			model.addConstr(f[k], GRB_LESS_EQUAL, f[k-1]);
		}
		
		// Set objective
		GRBLinExpr Obj = 0;
		for (int i = 1; i <= STP_N; ++i) { Obj += xIns.mComp[i].mCost * y[i]+STP_FIXEDCOST*u[i]; }
		model.setObjective(Obj, GRB_MINIMIZE);

		//initial solution
		vector<vector<double>> tempSol;
		double incumbent = ratio_heuristic(&xIns, tempSol); //has been checked
		for (int i = 1; i <= STP_N; i++) {
			for (int k = 1; k <= STP_N; k++) {
				x[i][k].set(GRB_DoubleAttr_Start, 0);
			}
		}
		for (int k = 0; k < tempSol.size(); k++) {
			for (int t = 4; t < tempSol[k].size(); t++) {
				x[int(tempSol[k][t])][k + 1].set(GRB_DoubleAttr_Start, 1);
			}
		}

		
		//Gurobi setting
		double time01 = double(double(clock() - startTime) / CLOCKS_PER_SEC);
		model.set(GRB_DoubleParam_TimeLimit, STP_SOVLERTIME - time01);
		//model.set(GRB_IntParam_Threads, 1);
		//model.set(GRB_IntParam_NumericFocus, 0);
		
		model.optimize(); 

		//cout <<fixed<< "Optimal: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
		
		
		// Get information
		//xIns.mObjVal = model.get(GRB_DoubleAttr_ObjVal);
		xIns.mState = model.get(GRB_IntAttr_Status);
		xIns.mGap = model.get(GRB_DoubleAttr_MIPGap);
		endTime = clock();
		if (xIns.mState > 2) { xIns.mCpuTime = STP_SOVLERTIME; }
		else { xIns.mCpuTime = double(double(endTime - startTime) / CLOCKS_PER_SEC); }
	
		compute_obj_val(xIns, x);
		
	}
	catch (GRBException ee) {
		cout << ee.getErrorCode() << endl;
		cout << ee.getMessage() << endl;

		xIns.mState = ee.getErrorCode();
		endTime = clock();
		xIns.mCpuTime = double(double(endTime - startTime) / CLOCKS_PER_SEC);
		cout << "This is Model P" << endl;
		cout << "Error Code: "<< ee.getErrorCode() << endl;
		system("pause");
	}
	catch (...) {
		cout << "Some problems happen" << endl;
		cout << "This is Model P" << endl;
		system("pause");
	}

		//Delete
		for (int i = 1; i <= STP_N; ++i) {delete[] x[i];}
		delete[] x; x = nullptr;

		delete[] y; y = nullptr;
		delete[] u; u = nullptr;
		delete[] f; f = nullptr;

		for (int i = 0; i <= STP_N; ++i) {delete[] z[i];}
		delete[] z; z = nullptr;

		return;
}




/*void assignment_based_model_validCut(StpInstance& xIns) {
	clock_t startTime, endTime;
	startTime = clock();
	//Define variable x_ik: (=1) if job i is assigned to batch k
	GRBVar** x = new GRBVar * [STP_N + 1];
	for (int i = 1; i <= STP_N; ++i) {
		x[i] = new GRBVar[STP_N + 1];
	}
	// Define variable y_i: the probability of component i to be executed 
	GRBVar* y = new GRBVar[STP_N + 1];

	//Define variable z_ik
	GRBVar** z = nullptr;
	z = new GRBVar * [STP_N + 1];
	for (int i = 0; i <= STP_N; ++i) {
		z[i] = new GRBVar[STP_N + 1];
	}

	// Define variable u_k: (=1) if the k-th batch is used
	GRBVar* u = new GRBVar[STP_N + 1];

	GRBVar* f = new GRBVar[STP_N + 1];
	try {
		GRBEnv env = GRBEnv();
		GRBModel model = GRBModel(env);

		//Define variable x_ik: (=1) if job i is assigned to batch k
		for (int i = 1; i <= STP_N; ++i) {
			for (int k = 1; k <= STP_N; ++k) {
				ostringstream vname;
				vname << "x_" << i << "_" << k;
				x[i][k] = model.addVar(0, 1, 0, GRB_BINARY, vname.str());
			}
		}
		// Define variable y_i: the probability of component i to be executed 
		for (int i = 1; i <= STP_N; ++i) {
			ostringstream vname;
			vname << "y_" << i;
			y[i] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, vname.str());
		}
		//Define variable z_ik
		for (int i = 0; i <= STP_N; ++i) {
			for (int k = 1; k <= STP_N; ++k) {
				ostringstream vname;
				vname << "z_" << i << "_" << k;
				z[i][k] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, vname.str());
			}
		}
		// Define variable u_k: (=1) if the k-th batch is used
		for (int k = 1; k <= STP_N; ++k) {
			ostringstream vname;
			vname << "u_" << k;
			u[k] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, vname.str());
		}

		for (int k = 1; k <= STP_N; ++k) {
			ostringstream vname;
			vname << "f_" << k;
			f[k] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, vname.str());
		}
		//////////////////////////////////////////
		//Add constraints
		//Component assignment
		for (int i = 1; i <= STP_N; i++) {
			GRBLinExpr expr = 0;
			for (int k = 1; k <= STP_N; ++k) {
				expr += x[i][k];
			}
			model.addConstr(expr, GRB_EQUAL, 1);
		}

		//Resource constraints
		for (int k = 1; k <= STP_N; k++) {
			GRBLinExpr expr = 0;
			for (int i = 1; i <= STP_N; i++) {
				expr += x[i][k] * xIns.mComp[i].mRes;
			}
			model.addConstr(expr, GRB_LESS_EQUAL, STP_RESCAP);
		}

		for (int i = 1; i <= STP_N; i++) {
			for (int k = 1; k <= STP_N; k++) {
				//GRBLinExpr expr = 0;
				//for (int v = 1; v <= k; v++) {
				//	expr += x[i][v];
				//}
				model.addConstr(y[i], GRB_GREATER_EQUAL, z[STP_N][k] - 1 + x[i][k]);
			}
		}

		for (int i = 0; i <= STP_N; i++) {
			model.addConstr(z[i][1], GRB_EQUAL, 1);
		}
		for (int k = 2; k <= STP_N; k++) {
			model.addConstr(z[0][k], GRB_EQUAL, z[STP_N][k - 1]);
		}

		for (int i = 1; i <= STP_N; i++) {
			for (int k = 2; k <= STP_N; k++) {
				model.addConstr(z[i][k], GRB_GREATER_EQUAL, z[i - 1][k] - x[i][k - 1]);
				model.addConstr(z[i][k], GRB_GREATER_EQUAL, xIns.mComp[i].mProb * z[i - 1][k]);
			}
		}

		for (int k = 1; k <= STP_N; k++) {
			for (int i = 1; i <= STP_N; i++) {
				model.addConstr(u[k], GRB_GREATER_EQUAL, z[STP_N][k] - 1 + x[i][k]);
			}
		}

		for (int k = 1; k <= STP_N; k++) {
			GRBLinExpr exp = 0;
			for (int i = 1; i <= STP_N; i++) {
				exp += x[i][k];
				model.addConstr(x[i][k], GRB_LESS_EQUAL, f[k]);
			}
			model.addConstr(exp, GRB_GREATER_EQUAL, f[k]);
		}
		for (int k = 2; k <= STP_N; k++) {
			model.addConstr(f[k], GRB_LESS_EQUAL, f[k - 1]);
		}

		for (int k = 2; k <= STP_N; k++){
			GRBLinExpr validCut = z[STP_N][k];
			for (int l = 1; l <= k - 1; ++l) {
				for (int i = 1; i <= STP_N; i++){
					validCut += x[i][l] * abs(log(xIns.mComp[i].mProb));
				}
			}
			model.addConstr(validCut, GRB_GREATER_EQUAL, 1);
		}

		// Set objective
		GRBLinExpr Obj = 0;
		for (int i = 1; i <= STP_N; ++i) { Obj += xIns.mComp[i].mCost * y[i] + STP_FIXEDCOST * u[i]; }

		//successive_knapsack(xIns, &env, &model, x);
		model.setObjective(Obj, GRB_MINIMIZE);
		double time01 = double(double(clock() - startTime) / CLOCKS_PER_SEC);
		model.set(GRB_DoubleParam_TimeLimit, STP_SOVLERTIME - time01);
		//model.set(GRB_IntParam_Threads, 1);


		model.optimize();

		// Get information
		//xIns.mObjVal = model.get(GRB_DoubleAttr_ObjVal);
		xIns.mState = model.get(GRB_IntAttr_Status);
		xIns.mGap = model.get(GRB_DoubleAttr_MIPGap);
		endTime = clock();
		if (xIns.mState > 2) { xIns.mCpuTime = STP_SOVLERTIME; }
		else { xIns.mCpuTime = double(double(endTime - startTime) / CLOCKS_PER_SEC); }

		compute_obj_val(xIns, x);

	}
	catch (GRBException ee) {
		cout << ee.getErrorCode() << endl;
		cout << ee.getMessage() << endl;

		xIns.mState = ee.getErrorCode();
		endTime = clock();
		xIns.mCpuTime = double(double(endTime - startTime) / CLOCKS_PER_SEC);
	}

	//Delete
	for (int i = 1; i <= STP_N; ++i) { delete[] x[i]; }
	delete[] x; x = nullptr;

	delete[] y; y = nullptr;
	delete[] u; u = nullptr;
	delete[] f; f = nullptr;

	for (int i = 0; i <= STP_N; ++i) { delete[] z[i]; }
	delete[] z; z = nullptr;

	return;
}
*/

void compute_obj_val(StpInstance &xIns, GRBVar** xVar) {
	vector<vector<int>> batchSeq;
	for (int k = 1; k <= STP_N; ++k) {
		vector<int> temp;
		for (int i = 1; i <= STP_N; i++) {
			if (xVar[i][k].get(GRB_DoubleAttr_X) > 0.5) {
				temp.push_back(i);
			}
		}
		if (temp.size() > 0) {batchSeq.push_back(temp);}
	}

	double totalCost = 0.0;
	double pp = 1.0;

	for ( int i = 0; i < batchSeq.size(); i++){
		double tempCost = double(STP_FIXEDCOST);
		for (int t = 0; t < batchSeq[i].size(); t++){
			tempCost += xIns.mComp[batchSeq[i][t]].mCost;
		}
		totalCost += tempCost * pp;

		for (int t = 0; t < batchSeq[i].size(); t++) {
			pp = pp * xIns.mComp[batchSeq[i][t]].mProb;
		}
	}
	xIns.mObjVal = totalCost;
	return;
}

void compute_obj_val(StpInstance& xIns, GRBVar** xVar,double &xObj) {
	vector<vector<int>> batchSeq;
	for (int k = 1; k <= STP_N; ++k) {
		vector<int> temp;
		for (int i = 1; i <= STP_N; i++) {
			if (xVar[i][k].get(GRB_DoubleAttr_X) > 0.5) {
				temp.push_back(i);
			}
		}
		if (temp.size() > 0) { batchSeq.push_back(temp); }
	}

	double totalCost = 0.0;
	double pp = 1.0;

	for (int i = 0; i < batchSeq.size(); i++) {
		double tempCost = STP_FIXEDCOST;
		for (int t = 0; t < batchSeq[i].size(); t++) {
			tempCost += xIns.mComp[batchSeq[i][t]].mCost;
		}
		totalCost += tempCost * pp;

		for (int t = 0; t < batchSeq[i].size(); t++) {
			pp = pp * xIns.mComp[batchSeq[i][t]].mProb;
		}
	}
	xObj = totalCost;
	return;
}






