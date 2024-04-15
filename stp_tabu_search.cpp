#include "stp_tabu_search.h"
#include "stp_ratio_heuristic.h"
#include "stp_model_solver.h"

int createArcs(StpInstance *xIns, vector<vector<int> >& arcs, vector<bool>& souvenir) {

	int m = STP_N;      // component number
	int c = STP_RESCAP; // resource constraint
	vector<int> w; w.push_back(-1);      // component resource requirement
	vector <int> b; b.push_back(-1);     // component 

	for (int i = 1; i <= STP_N; i++){
		w.push_back(xIns->mComp[i].mRes);
		b.push_back(1);
	}

	// Local variables
	souvenir.resize(c + 1, false);
	souvenir[0] = true;
	vector<int> vec(3);
	set<vector<int> > setArcs;

	// Arcs for items		
	for (int i = 1; i <= m; i++) {
		for (int j = c - w[i]; j >= 0; j--) {
			if (souvenir[j] == true) {
				for (int k = 1; k <= b[i]; k++) {
					if (j + k * w[i] > c) {
						break;
					}
					souvenir[j + k * w[i]] = true;
					vec[0] = j + (k - 1) * w[i];
					vec[1] = j + k * w[i];
					vec[2] = i;
					setArcs.insert(vec);
				}
			}
		}
	}

	// Loss Arcs
	vector<int> activeN;
	for (int i = 1; i <= c; i++) {
		if (souvenir[i] == true) {
			activeN.push_back(i);
		}
	}
	for (int i = 0; i < activeN.size() - 1; i++) {
		vec[0] = activeN[i]; vec[1] = activeN[i + 1]; vec[2] = -1;
		setArcs.insert(vec);
	}
	//add the last loss arc, when there is no arc to connect the terminal node
	if (!activeN.empty()) {
		vec[0] = activeN.back(); vec[1] = c; vec[2] = -1;
		setArcs.insert(vec);
	}

	for (set<vector<int> >::iterator it = setArcs.begin(); it != setArcs.end(); ++it) {
		arcs.push_back(*it);
	}

	return 0;
}

int transformSolution(vector<vector<vector<int> > >& arcsUsed, vector<vector<int> >& sol) {
	int c = STP_RESCAP;
	sol.resize(0);
	while (arcsUsed[0].size() > 0) {
		vector<int> curBin;
		int thisTail = 0;
		for (;;) {

			int nextTail = arcsUsed[thisTail].back()[1];

			// If the next arc is an item
			if (arcsUsed[thisTail].back()[2] >= 0) {
				curBin.push_back(arcsUsed[thisTail].back()[2]);
			}

			arcsUsed[thisTail].pop_back();
			thisTail = nextTail;

			if (thisTail == c) {
				break;
			}
		}
		sol.push_back(curBin);
	}
	return 0;
}

int arcflow(StpInstance *xIns, vector<vector<int> >& sol) {
	clock_t startTime;
	startTime = clock();
	
	vector<double> info; info.resize(7);
	int of = 9999999;  //Objective value
	int m = STP_N;      // component number
	int c = STP_RESCAP; // resource constraint
	vector<int> w; w.push_back(-1);      // component resource requirement
	vector <int> b; b.push_back(-1);     // component 

	for (int i = 1; i <= STP_N; i++) {
		w.push_back(xIns->mComp[i].mRes);
		b.push_back(1);
	}
	
	
	// Local variables
	vector<vector<int> > arcs;
	vector<bool> souvenir;

	// Arc definition
	createArcs(xIns, arcs, souvenir);

	// Fill info
	info[1] = double(double(clock() - startTime) / CLOCKS_PER_SEC);
	info[2] = double(arcs.size());
	info[3] = 0;
	for (int i = 0; i < c; i++) {
		if (souvenir[i] == true) {
			info[3]++;
		}
	}

	if (arcs.size() > 2000000) {
		info[0] = double(double(clock() - startTime) / CLOCKS_PER_SEC);
		return -1;
	}

	// Model
	GRBEnv env = GRBEnv();
	try {
		// Close the output about setting parameter
		env.set(GRB_IntParam_OutputFlag, 0);
		GRBModel model = GRBModel(env);

		// Declaration of the variables for the model
		vector<GRBVar> isArcUsed(arcs.size());
		vector<GRBLinExpr>	cIn(c + 1);
		vector<GRBLinExpr>  cOut(c + 1);
		vector<GRBLinExpr>  objT(m+1);

		// Initialization
		for (int i = 0; i < arcs.size(); i++) {
			if (arcs[i][2] >= 0) {
				isArcUsed[i] = model.addVar(0, GRB_INFINITY, 0, GRB_INTEGER);
			}
			else {
				isArcUsed[i] = model.addVar(0, GRB_INFINITY, 0, GRB_INTEGER);
			}
		}
		for (int i = 0; i < c + 1; i++) {
			cIn[i] = 0;
			cOut[i] = 0;
		}
		for (int i = 1; i <= m; i++) {
			objT[i] = 0;
		}

		model.update();

		// Constraints
		for (int i = 0; i < arcs.size(); i++) {
			cIn[arcs[i][1]] += isArcUsed[i];
			cOut[arcs[i][0]] += isArcUsed[i];
			if (arcs[i][2] >= 0) {
				objT[arcs[i][2]] += isArcUsed[i];
			}
		}

		// Objective function
		model.setObjective(cOut[0], GRB_MINIMIZE);
		
		// Flow conservation
		for (int i = 1; i < c; i++) {
			if (souvenir[i] == true) {
				model.addConstr(cIn[i] == cOut[i]);
			}
		}
		
		model.addConstr(cOut[0] == cIn[c]);
		
		// Take the objects
		for (int i = 1; i <= m; i++) {
			if (b[i] == 1) {
				model.addConstr(objT[i] == b[i]);
			}
			else {
				model.addConstr(objT[i] == b[i]);
			}
		}
		
		// Setting of Gurobi
		//double time01 = double(double(clock() - startTime) / CLOCKS_PER_SEC);
		//model.set(GRB_DoubleParam_TimeLimit, STP_SOVLERTIME - time01);
		//model.set(GRB_IntParam_Threads, 1);
		model.set(GRB_IntParam_LogToConsole, 0);
	
		model.optimize();

		// Fill info
		info[4] = model.get(GRB_IntAttr_NumVars);
		info[5] = model.get(GRB_IntAttr_NumConstrs);
		info[6] = model.get(GRB_IntAttr_NumNZs);

		// Solving the model (or get some information)
		if (model.get(GRB_IntAttr_SolCount) < 1) {
			cout << "Failed to optimize LP. " << endl;
			// get info
			info[0] = double(double(clock() - startTime) / CLOCKS_PER_SEC);
			return -1;
		}

		// Get solution
		vector<vector<vector<int> > > arcsUsed(c + 1);
		for (int i = 0; i < arcs.size(); i++) {
			//if(isArcUsed[i].get(GRB_DoubleAttr_X)>0.5 )
			//cout << isArcUsed[i].get(GRB_DoubleAttr_X)<<"  "<<arcs[i][2] << endl;
			for (int j = 0; j < ceil(isArcUsed[i].get(GRB_DoubleAttr_X)); j++) {
				arcsUsed[arcs[i][0]].push_back(arcs[i]);
			}
		}

		// Fill info
		info[0] = double(double(clock() - startTime) / CLOCKS_PER_SEC);
		if (of >= ceil(model.get(GRB_DoubleAttr_ObjVal))) {
			of = ceil(model.get(GRB_DoubleAttr_ObjVal));
			//xIns->mObjVal = of;
			transformSolution(arcsUsed, sol);
		}
	}

	// Exceptions
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		info[0] = double(double(clock() - startTime) / CLOCKS_PER_SEC);
		cout << e.getMessage() << endl;
		system("pause");
	}
	catch (...) {
		cout << "Exception during optimization" << endl;
		info[0] = double(double(clock() - startTime) / CLOCKS_PER_SEC);
		system("pause");
	}

	// End
	return 0;
}

void binPacking(StpInstance* xIns) {
	clock_t startTime, endTime;
	startTime = clock();
	
	//Define variable x_ik: (=1) if job i is assigned to batch k
	GRBVar** x = new GRBVar * [STP_N + 1];
	for (int i = 1; i <= STP_N; ++i) {
		x[i] = new GRBVar[STP_N + 1];
	}
	
	// Define variable u_k: (=1) if the k-th batch is used
	GRBVar* u = new GRBVar[STP_N + 1];

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

		// Define variable u_k: (=1) if the k-th batch is used
		for (int k = 1; k <= STP_N; ++k) {
			ostringstream vname;
			vname << "u_" << k;
			u[k] = model.addVar(0, GRB_INFINITY, 0, GRB_BINARY, vname.str());
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
				expr += x[i][k] * xIns->mComp[i].mRes;
			}
			model.addConstr(expr, GRB_LESS_EQUAL, STP_RESCAP*u[k]);
		}

		// Set objective
		GRBLinExpr Obj = 0;
		for (int k = 1; k <= STP_N; ++k) { Obj += u[k]; }

	
		model.setObjective(Obj, GRB_MINIMIZE);
		double time01 = double(double(clock() - startTime) / CLOCKS_PER_SEC);
		model.set(GRB_DoubleParam_TimeLimit, STP_SOVLERTIME - time01);
		//model.set(GRB_IntParam_Threads, 1);


		model.optimize();

		// Get information
		xIns->mObjVal = model.get(GRB_DoubleAttr_ObjVal);
		xIns->mState = model.get(GRB_IntAttr_Status);
		endTime = clock();
		xIns->mCpuTime = double(double(endTime - startTime) / CLOCKS_PER_SEC);
		xIns->mGap = model.get(GRB_DoubleAttr_MIPGap);

	}
	catch (GRBException ee) {
		cout << ee.getErrorCode() << endl;
		cout << ee.getMessage() << endl;

		xIns->mState = ee.getErrorCode();
		endTime = clock();
		xIns->mCpuTime = double(double(endTime - startTime) / CLOCKS_PER_SEC);
	}

	//Delete
	for (int i = 1; i <= STP_N; ++i) { delete[] x[i]; }
	delete[] x; x = nullptr;
	delete[] u; u = nullptr;

	return;
}

TabuSolver::TabuSolver(StpInstance *xIns) {
	startTime = clock();
	ins = xIns;
	counterB = 0;

	//initial_solution(ins);
	memset(tabuList, 0, sizeof(int)*(STP_N+1));
	memset(visit, 0, sizeof(int) * (STP_N + 1));
	visit[0] = STP_VMAX + 1;
	currIteration = 0;
	compMoved = 0;
	repairMoved = 0;
	tabuStop = 100;
	//for re-initialization
	reIter = 0;
	
	final_bestObj = double(INT_MAX);
	final_bestSol.resize(0);
	//std::random_device rd;
	//std::mt19937 gen(rd());
	//std::uniform_int_distribution<> dis(5, 15);
	//tabuTenure =  dis(gen);
	
	//Initialize best obj & solution
	/*for (int i = 1; i <= STP_N; i++) {
		vector<double> tempBatch = { double(ins->mComp[i].mRes),0.0,
			ins->mComp[i].mProb,double(STP_FIXEDCOST+ ins->mComp[i].mCost),double(i) };
		tempBatch[1] = tempBatch[3] / (1 - tempBatch[2]);
		this->bestSol.push_back(tempBatch);
	}
	this->compute_ObjVal(bestObj, bestSol);
	*/

	// For hybrid proximity search
	mMode = 0;
	mVar = nullptr;
	mDelta = nullptr;
	mIncumbent = nullptr;
	mSol = nullptr;
	mImproved = nullptr;
	
	return;
}

void TabuSolver::initial_solution(StpInstance* xIns) {

	vector<vector<int>> sol;
	arcflow(xIns, sol);
	this->batchLB = int(sol.size());
	for (int i = 0; i < sol.size(); i++) {
		vector<double> tempBatch = { 0.0,0.0,1.0,STP_FIXEDCOST };
		for (int v = 0; v < sol[i].size(); v++) {
			tempBatch.push_back(sol[i][v]);
			tempBatch[0] += xIns->mComp[sol[i][v]].mRes;
			tempBatch[2] = tempBatch[2] * xIns->mComp[sol[i][v]].mProb;
			tempBatch[3] += xIns->mComp[sol[i][v]].mCost;
		} 
		tempBatch[1] = tempBatch[3] / (1 - tempBatch[2]);
		this->currSol.push_back(tempBatch);
	}

	compute_ObjVal(this->currObj,currSol);
};

bool LessSort(vector<double>& a, vector<double>& b) { return (a[1] < b[1]); }

/*void TabuSolver::compute_ObjVal(double& xObj, vector<vector<double>>& xBatchSeq) {
	// Caculate objective value
	sort(xBatchSeq.begin(), xBatchSeq.end(), LessSort);
	double prob= 1.0;
	xObj = 0.0;
	for (int i = 0; i < xBatchSeq.size(); i++) {
		xObj += xBatchSeq[i][3] * prob;
		prob = prob * xBatchSeq[i][2];
	}
	
	return;
}*/

void TabuSolver::get_best_neighbor() {

	vector<vector<double>> currCopy = currSol;
	double objCopy = currObj;

	vector<vector<double>> currNeighbor = currSol;
	double objNeighbor = INT_MAX;
	
	for (int i = 0; i < currCopy.size(); i++){
		//if (currCopy[i].size() > 5) {
			for (int t = 4; t < currCopy[i].size(); t++) {
				//if non-tabu
				if (tabuList[int(currCopy[i][t])] < currIteration) {
					for (int k = 0; k < currCopy.size(); k++) {
						if (k != i) {
							vector<vector<double>> tempSol = currSol;
							double tempObj = INT_MAX;
							if (currCopy[k][0] + ins->mComp[int(currCopy[i][t])].mRes <= STP_RESCAP) {

								// insertion the component to another batch
								comp_insertion(int(currCopy[i][t]), i, k, tempSol);
								compute_ObjVal(tempObj, tempSol);
							
								//if (tempObj < objNeighbor) {
								if (objNeighbor-tempObj>STP_EPSILON) {
									objNeighbor = tempObj;
									currNeighbor = tempSol;
									compMoved = int(currCopy[i][t]);
									repairMoved = 0;
									ins->tempTTB03 = double(double(clock() - startTime) / CLOCKS_PER_SEC);
								}
								//cout << "tempObj: " << tempObj << " Neighbor: " << objNeighbor << endl;

							}
							else {
								//repair
								int deltaR = int(currCopy[k][0]) + ins->mComp[int(currCopy[i][t])].mRes - STP_RESCAP;
								int tempRepairComp = 0;
								bool repaired = false;
								for (int v = 4; v < currCopy[k].size(); v++) {
									if (ins->mComp[int(currCopy[k][v])].mRes >= deltaR
										&& (tabuList[int(currCopy[k][v])] < currIteration)
										) {
										// find a batch to move the component u
										for (int g = 0; g < currCopy.size(); g++) {
											if ((g != k) && (g != i)
												&& (currCopy[g][0] + ins->mComp[int(currCopy[k][v])].mRes <= STP_RESCAP))
											{
												//cout << "repaired by: "<< currCopy[k][v] << endl;
												// for other batch, which is not i or k
												//tabuList[int(currCopy[k][v])] = currIteration + tabuTenure;
												comp_insertion(int(currCopy[k][v]), k, g, tempSol);
												tempRepairComp = int(currCopy[k][v]);
												repaired = true;
												break;
											}
											else if ((g == i)
												&& (currCopy[g][0] + ins->mComp[int(currCopy[k][v])].mRes - ins->mComp[int(currCopy[i][t])].mRes <= STP_RESCAP)
												) // for batch i, similar to swap
											{
												//cout << "repaired by: " << currCopy[k][v] << endl;
												//tabuList[int(currCopy[k][v])] = currIteration + tabuTenure;
												comp_insertion(int(currCopy[k][v]), k, g, tempSol);
												tempRepairComp = int(currCopy[k][v]);
												repaired = true;
												break;
											}
										}
										if (repaired) { break; }
									}
								}
								//if repaired, then move the component
								if (repaired) {

									comp_insertion(int(currCopy[i][t]), i, k, tempSol);
									compute_ObjVal(tempObj, tempSol);
									//if (tempObj < objNeighbor) {
									if (objNeighbor - tempObj > STP_EPSILON) {
										objNeighbor = tempObj;
										currNeighbor = tempSol;
										compMoved = int(currCopy[i][t]);
										repairMoved = tempRepairComp;
										ins->tempTTB03 = double(double(clock() - startTime) / CLOCKS_PER_SEC);
									}
								}
								//else {cout << "No feasible batch to repair" << endl;}
							}

							

						}
					}

				}
			}
		//}
	}

	//cout << "last: " << currObj << "  Neighbor: " << objNeighbor << endl;
	if (compMoved > 0) {
		currSol = currNeighbor;
		currObj = objNeighbor;
	}
	//cout << "Moved: " << compMoved<<" repaired: "<<repairMoved << endl;
	return;
}

void TabuSolver::comp_insertion(int xComp,  int xOriBatch,
	int xDestBatch, vector<vector<double>>& xBatchSeq) {

	//for destination batch
	xBatchSeq[xDestBatch][0] += ins->mComp[xComp].mRes;
	xBatchSeq[xDestBatch][3] += ins->mComp[xComp].mCost;
	xBatchSeq[xDestBatch][2] = xBatchSeq[xDestBatch][2] * ins->mComp[xComp].mProb;
	xBatchSeq[xDestBatch][1] = xBatchSeq[xDestBatch][3] / (1 - xBatchSeq[xDestBatch][2]);
	xBatchSeq[xDestBatch].push_back(double(xComp));

	//for original batch
	vector<double>::iterator it = find(xBatchSeq[xOriBatch].begin() + 4, xBatchSeq[xOriBatch].end(), double(xComp));
	if (it != xBatchSeq[xOriBatch].end()){
		xBatchSeq[xOriBatch].erase(it);
	}
	else {
		cout << "No corresponding component" << endl;
		system("pause");
	}
	xBatchSeq[xOriBatch][0] -= ins->mComp[xComp].mRes;
	xBatchSeq[xOriBatch][3] -= ins->mComp[xComp].mCost;
	if (xBatchSeq[xOriBatch].size() > 4) {
		xBatchSeq[xOriBatch][2] = xBatchSeq[xOriBatch][2] / ins->mComp[xComp].mProb;
		xBatchSeq[xOriBatch][1] = xBatchSeq[xOriBatch][3] / (1 - xBatchSeq[xOriBatch][2]);
	}
	else {
		xBatchSeq[xOriBatch][2] = 1;
		xBatchSeq[xOriBatch][1] = INT_MAX;
		xBatchSeq[xOriBatch][0] = 0.0;
		xBatchSeq[xOriBatch][3] = STP_FIXEDCOST;
	}
	
	return;
}

void TabuSolver::basic_tabu_search() {

	memset(visit, 0, sizeof(int) * (STP_N + 1));
	visit[0] = STP_VMAX + 1;
	while (counterB <= STP_N && double(double(clock() - startTime) / CLOCKS_PER_SEC) < STP_SOVLERTIME)
	{
		compute_ObjVal(this->currObj, currSol);
		//cout << "batch LB: " << batchLB << " this iteration:  " << counterB << endl;
		visit[counterB]++;
		
		

		memset(tabuList, 0, sizeof(int) * (STP_N + 1));
		currIteration = 1;
		int stopCriteria = 0;

		while (stopCriteria < this->tabuStop
			&& double(double(clock() - startTime) / CLOCKS_PER_SEC) < STP_SOVLERTIME)
		{
			//generate a new non-tabu solution
			compMoved = 0;
			repairMoved = 0;
			if (double(double(clock() - startTime) / CLOCKS_PER_SEC) < STP_SOVLERTIME)
			{
			
				get_best_neighbor();
			}
			else { break; }
			// update tabu list
			if (compMoved > 0) { tabuList[compMoved] = currIteration + tabuTenure; }
			if (repairMoved > 0) { tabuList[repairMoved] = currIteration + tabuTenure; }
			

			if (bestObj-currObj>STP_EPSILON){
				ins->tempTTB02 = ins->tempTTB03;
				//cout <<fixed<< "before: " << bestObj << "  after: " << currObj << endl;
				bestSol = currSol; bestObj = currObj;
				stopCriteria = 0;
				//cout << "Update!" << endl;
			}
			else { stopCriteria++; }
			currIteration++;
		}
		
		if (!batch_reduction()) {
			do {
				// Add an empty batch
				vector<double> temp = { 0,double(INT_MAX),1.0,STP_FIXEDCOST };
				currSol.push_back(temp);
				compute_ObjVal(this->currObj, currSol);
				counterB++;
			} while (visit[counterB] >= STP_VMAX && counterB <= STP_N);
			//cout << "batch reduction fails, Now batch number: " << counterB << endl;
		}
		else {counterB--;
			//cout << "Batch reduction success: " << counterB << endl; 
		}
		if (bestObj - currObj > STP_EPSILON) { 
			ins->tempTTB02 = double(double(clock() - startTime) / CLOCKS_PER_SEC);
			bestSol = currSol; bestObj = currObj; }
		//cout << endl;
		/*// Add an empty batch
		vector<double> temp = { 0,INT_MAX,1,STP_FIXEDCOST };
		currSol.push_back(temp);
		this->compute_ObjVal(this->currObj, currSol);
		counterB++;*/
	}

	return;
}

bool TabuSolver::batch_reduction() {
	if (visit[int(currSol.size()-1)]>= STP_VMAX 
		|| (int(currSol.size())== batchLB)) {
		
		return false;
	}
	vector<int> rMax(currSol.size(), -1);
	// find out the max resource requirement of each batch
	for (int i = 0; i < currSol.size(); i++){
		for (int v = 4; v < currSol[i].size(); v++){
			if (ins->mComp[int(currSol[i][v])].mRes > rMax[i]) {
				rMax[i] = ins->mComp[int(currSol[i][v])].mRes;
			}
		}
		
	}
	
	//int pos = int(std::min_element(rMax.begin(), rMax.end())-rMax.begin());
	int pos = weighed_random(rMax);
	vector<int>  list;
	for (int i = 4; i < currSol[pos].size(); i++){
		list.push_back(int(currSol[pos][i]));
	}
	//cout << "The batch to be removed has job: " <<list.size()<< endl;

	// order job resouce by a decreasing order
	for (int i = 0; i < list.size(); i++){
		int min=i;
		for (int j = i+1; j <list.size(); j++){
			if (ins->mComp[list[j]].mRes > ins->mComp[list[min]].mRes) {
				min = j;
			}
		}
		int temp = list[min];
		list[min] = list[i];
		list[i] = temp;
	}

	vector<vector<double>> copySol = currSol;

	for (int i = 0; i < list.size(); i++){
		int destination = -1;
		//int remainRes = INT_MAX;
		for (int k = 0; k < currSol.size(); k++){
			if (k != pos && currSol[k][0]+ins->mComp[list[i]].mRes<=STP_RESCAP)
				//&& STP_RESCAP- currSol[k][0]- ins->mComp[list[i]].mRes<remainRes) 
			{
				destination = k;
				break;
				//remainRes = int(STP_RESCAP - currSol[k][0] - ins->mComp[list[i]].mRes);
			}
		}
		if (destination < 0) {
			currSol = copySol;
			return false;}
		else {comp_insertion(list[i], pos, destination, currSol);}
	}

	
	// Order solution & delete empty batch
	sort(currSol.begin(), currSol.end(), LessSort);
	currSol.pop_back();
	compute_ObjVal(currObj, currSol);
	

	return true;
}

int weighed_random(vector<int>& xVec) {
	
	double sum = 0;
	vector<double> weight;
	vector<double> v;

	for (int i = 0; i < xVec.size(); i++) {
		if (xVec[i] > 0) {
			weight.push_back(1.0 / xVec[i]);
			sum += weight[i];
		}
	}

	for (int i = 0; i < weight.size(); i++) {
		double temp = 0.0;
		for (int j = 0; j <= i; j++) {
			temp += weight[j];
		}
		v.push_back(temp);
	}

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0, sum);
	double randomNumber = dis(gen); 

	for (int i = 0; i < v.size(); ++i) {
		if (randomNumber <= v[i])
			return i;
	}
};

TabuSolver::TabuSolver(StpInstance* xIns,GRBVar** xVar,
	GRBLinExpr* xDelta, double* xIncumbent, int xMode, double** xSol, bool* xImproved, clock_t xStart) {
	// For Hybrid Proximity
	mMode = xMode;
	mVar = xVar;
	mDelta = xDelta;
	mIncumbent = xIncumbent;
	mSol = xSol;
	mImproved = xImproved;

	counterB = 0;
	startTime = xStart;
	ins = xIns;
	if (mMode < 2) { 
		tabuStop = 50; 
		final_bestObj = double(INT_MAX);
		final_bestSol.resize(0);
	}
	else {
		tabuStop = 50;
		// compute currSol & currObj
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
		currObj = 0.0;
		currSol.clear();
		//this->batchLB = int(batchSeq.size());
		for (int i = 0; i < batchSeq.size(); i++) {
			vector<double> tempBatch = { 0.0,0.0,1.0,STP_FIXEDCOST };
			for (int v = 0; v < batchSeq[i].size(); v++) {
				tempBatch.push_back(batchSeq[i][v]);
				tempBatch[0] += xIns->mComp[batchSeq[i][v]].mRes;
				tempBatch[2] = tempBatch[2] * xIns->mComp[batchSeq[i][v]].mProb;
				tempBatch[3] += xIns->mComp[batchSeq[i][v]].mCost;
			}
			tempBatch[1] = tempBatch[3] / (1 - tempBatch[2]);
			this->currSol.push_back(tempBatch);
		}

		compute_ObjVal(this->currObj, currSol);
		final_bestObj = currObj;
		final_bestSol = currSol;
	}

	memset(tabuList, 0, sizeof(int) * (STP_N + 1));
	memset(visit, 0, sizeof(int) * (STP_N + 1));
	visit[0] = STP_VMAX + 1;

	currIteration = 0;
	compMoved = 0;
	repairMoved = 0;
	reIter = 0;
	
	
	
	return;

}



void TabuSolver::solve_re_initial() {
	cout << "Tabu Search starts" << endl;
	counterB = 0;
	while (double(double(clock() - startTime) / CLOCKS_PER_SEC) < STP_SOVLERTIME) {
		// differen intialization
		if (reIter == 0 && this->mMode <= 1) {
			cout << "reIter: " << reIter << endl;
			this->initial_solution(ins);
			counterB = int(currSol.size());
		}
		else if (reIter == 0 && this->mMode > 1) {
			counterB = int(currSol.size());
		}
		else if (reIter == 1) {
			ratio_heuristic(ins, this->currSol); //has been checked
			counterB = int(this->currSol.size());
			compute_ObjVal(this->currObj, currSol);
			cout << "reIter: " << reIter << endl;
		}
		else {
			//cout << "reIter: " << reIter << endl;
			random_initialization(ins, this->currSol);
			counterB = int(this->currSol.size());
			compute_ObjVal(this->currObj, currSol);
		}

		bestObj = currObj;
		bestSol = currSol;

		if (final_bestObj - bestObj > STP_EPSILON) {
			ins->tempTTB = double(double(clock() - startTime) / CLOCKS_PER_SEC);
			final_bestObj = bestObj;
			final_bestSol = bestSol;
		}

		//initialize tabu tenure
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_int_distribution<> dis(5, 15);
		tabuTenure = dis(gen);

		


		this->basic_tabu_search();

		if (final_bestObj - bestObj > STP_EPSILON) {
			if (mMode < 2) { ins->tempTTB = ins->tempTTB02; }
			final_bestObj = bestObj;
			final_bestSol = bestSol;
		}

		reIter++;
		if (reIter > 1 && this->mMode == 1) {
			cout << "break from mode 1, re_iter: " << reIter << endl;
			break;
		}
		else if (reIter > 0 && this->mMode == 2) {
			cout << "break from mode 2, re_iter: " << reIter << endl;
			break;
		}

	}

	if (this->mMode == 1) {
		// set initial solution for hybrid proximity search
		if (*mIncumbent - final_bestObj > STP_EPSILON) {
			*this->mDelta = 0;

			for (int i = 1; i <= STP_N; i++) {
				for (int k = 1; k <= STP_N; k++) {
					this->mVar[i][k].set(GRB_DoubleAttr_Start, 0);
					*this->mDelta += mVar[i][k];
				}
			}
			for (int k = 0; k < final_bestSol.size(); k++) {
				for (int t = 4; t < final_bestSol[k].size(); t++) {
					this->mVar[int(final_bestSol[k][t])][k + 1].set(GRB_DoubleAttr_Start, 1);
					*this->mDelta += (1 - this->mVar[int(final_bestSol[k][t])][k + 1] - this->mVar[int(final_bestSol[k][t])][k + 1]);
				}
			}
			//mEta->set(GRB_DoubleAttr_Start, final_bestObj);
			*mIncumbent = final_bestObj;
		}
		return;
	}
	else if (mMode == 2) {
		if (*mIncumbent - final_bestObj > STP_EPSILON) {
			*mIncumbent = final_bestObj;
			*mImproved = true;
			ins->tempTTB = ins->tempTTB02;
			for (int k = 0; k < final_bestSol.size(); k++) {
				for (int t = 4; t < final_bestSol[k].size(); t++) {
					mSol[int(final_bestSol[k][t])][k + 1] = 1;
				}
			}
		}
		return;
	}

	this->endTime = clock();
	ins->mObjVal = final_bestObj;
	ins->mCpuTime = double(double(endTime - startTime) / CLOCKS_PER_SEC);
	ins->timeToBest = ins->tempTTB;
	return;
}

void hybrid_proximity(StpInstance& xIns) {
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
	double incumbent = double(INT_MAX);

	try {
		GRBEnv env = GRBEnv();
		GRBModel model = GRBModel(env);
		GRBLinExpr delta = 0; //record Hamming distance between solutions
		double objValue = 0.0;
		
		GRBVar flexVar = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
		double** xSol = new double* [STP_N + 1];
		for (int i = 1; i <= STP_N; i++) {xSol[i] = new double[STP_N + 1] { 0 };}
		bool improved = false;
		
		TabuSolver tabu(&xIns, x,  &delta, &incumbent, 1, xSol, &improved, startTime);

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

		// Obj function expression
		GRBLinExpr objFunction = 0;
		for (int i = 1; i <= STP_N; ++i) { objFunction += xIns.mComp[i].mCost * y[i] + STP_FIXEDCOST * u[i]; }

		// Set parameters
		double step = 0.01;
		double bigM =  10000;
		

		// initial reference solution
		tabu.solve_re_initial();
		cout << endl;  
		cout << "Initial objective value by Tabu: " << incumbent << endl; 
		cout << endl;


		//Gurobi Setting
		//model.set(GRB_IntParam_Threads, 1);
		model.set(GRB_DoubleParam_MIPGapAbs, bigM * 0.2);
		//model.set(GRB_IntParam_LogToConsole,0);
		int counter = 0;
		int stop = 0;
		while (((double(clock() - startTime) / CLOCKS_PER_SEC) < STP_SOVLERTIME) && stop < 10) {

			double time01 = STP_SOVLERTIME - double(double(clock() - startTime) / CLOCKS_PER_SEC);
			model.set(GRB_DoubleParam_TimeLimit, min_number(double(STP_N), time01));

			//Set cut off constraint and objective
			model.addConstr(objFunction <= incumbent - step * incumbent * (1 - flexVar), "Cut_off");
			model.setObjective(delta + flexVar * bigM, GRB_MINIMIZE);
			//model.set(GRB_IntParam_NumericFocus,2);
			//model.set(GRB_IntParam_SolutionLimit, 3);

			model.optimize();

			objValue = 0.0;
			//compute_obj_val(xIns, x, objValue);
			//cout << endl;
			//cout << "++++  Refine" << endl;
			refine(xIns, x, &objValue, startTime);
			cout <<fixed<< " Objective value by proximity search:  " << objValue << endl;

			//Tabu parameter setting
			for (int i = 1; i <= STP_N; i++) { memset(xSol[i], 0, sizeof(double) * (STP_N + 1)); }
			improved = false;


			if (STP_EPSILON < incumbent - objValue) {
				cout << "Proximity search improved " << endl;
				cout << "Incumbent: " << incumbent << "  Proximity Search: " << objValue << "  update " << endl;
				
				incumbent = objValue;
				xIns.tempTTB = double(double(clock() - startTime) / CLOCKS_PER_SEC);

				TabuSolver tabu02(&xIns, x, &delta, &incumbent, 2, xSol, &improved, startTime);
				tabu02.batchLB = tabu.batchLB;
				tabu02.solve_re_initial();
				
				cout << "Tabu Improved: " << improved << endl;
				cout << "incumbent after Tabu: " << incumbent << endl;

				model.remove(model.getConstrByName("Cut_off"));
				//eta.set(GRB_DoubleAttr_Start, incumbent);
				delta = 0;
				if (improved) {	
					for (int i = 1; i <= STP_N; ++i) {
						for (int k = 1; k <= STP_N; ++k) {
							if (xSol[i][k] > 0.5) {
								delta += 1 - x[i][k];
								x[i][k].set(GRB_DoubleAttr_Start, 1);
							}
							else {
								delta += x[i][k];
								x[i][k].set(GRB_DoubleAttr_Start, 0);
							}
						}
					}
					
				}
				else {
					for (int i = 1; i <= STP_N; ++i) {
						for (int k = 1; k <= STP_N; ++k) {
							if (x[i][k].get(GRB_DoubleAttr_X) > 0.5) {
								delta += 1 - x[i][k];
								x[i][k].set(GRB_DoubleAttr_Start, 1);
							}
							else {
								delta += x[i][k];
								x[i][k].set(GRB_DoubleAttr_Start, 0);
							}
						}
					}
				}

			}
			else {
				cout << endl;
				cout << "adjust step" << endl;
				model.remove(model.getConstrByName("Cut_off"));
				step = step * 0.8;
				cout<< "step: " << step*incumbent << endl;
				stop++;
				TabuSolver tabu02(&xIns, x,  &delta, &incumbent, 2, xSol, &improved, startTime);
				tabu02.batchLB = tabu.batchLB;
				tabu02.solve_re_initial();
				
				cout << "Improved after Tabu: " << improved << endl;
				cout << "Incumbent after Tabu: " << incumbent << endl;
				cout << endl;

				if (improved) {
					//eta.set(GRB_DoubleAttr_Start, incumbent);
					delta = 0;
					for (int i = 1; i <= STP_N; ++i) {
						for (int k = 1; k <= STP_N; ++k) {
							if (xSol[i][k] > 0.5) {
								delta += 1 - x[i][k];
								x[i][k].set(GRB_DoubleAttr_Start, 1);
							}
							else {
								delta += x[i][k];
								x[i][k].set(GRB_DoubleAttr_Start, 0);
							}
						}
					}
				}
			}
			counter++;
			cout << "Iteration: " << counter << " No improvements: " << stop << endl;
			cout << endl;
		}


		// Get information
		xIns.mObjVal = incumbent;
		endTime = clock();
		xIns.mCpuTime = double(double(endTime - startTime) / CLOCKS_PER_SEC);
		//xIns.mGap = model.get(GRB_DoubleAttr_MIPGap);
		xIns.timeToBest = xIns.tempTTB;

		for (int i = 1; i <= STP_N; i++) { delete[] xSol[i]; }
		delete[] xSol;

	}
	catch (GRBException ee) {
		cout << ee.getErrorCode() << endl;
		cout << ee.getMessage() << endl;
		xIns.mObjVal = incumbent;
		xIns.mState = ee.getErrorCode();
		endTime = clock();
		xIns.mCpuTime = double(double(endTime - startTime) / CLOCKS_PER_SEC);
	}
	catch (...) {
		cout << "Problems Happen. This is Hybrid Proximity" << endl;
		system("pause");
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
};


void refine(StpInstance& xIns, GRBVar** xVar, double* xObjValue,clock_t xStart) {
	clock_t startTime, endTime;
	startTime = xStart;
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
				if (xVar[i][k].get(GRB_DoubleAttr_X) > 0.5){model.addConstr(x[i][k], GRB_EQUAL, 1);}
				else { model.addConstr(x[i][k], GRB_EQUAL, 0); }
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
				//for (int v = 1; v <=k; v++){
					//expr += x[i][v];
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

		// Set objective
		GRBLinExpr Obj = 0;
		for (int i = 1; i <= STP_N; ++i) { Obj += xIns.mComp[i].mCost * y[i] + STP_FIXEDCOST * u[i]; }
		model.setObjective(Obj, GRB_MINIMIZE);


		//Gurobi setting
		//double time01 = double(double(clock() - startTime) / CLOCKS_PER_SEC);
		//model.set(GRB_DoubleParam_TimeLimit, STP_SOVLERTIME - time01);
		//model.set(GRB_IntParam_Threads, 1);
		//model.set(GRB_IntParam_NumericFocus, 0);
		model.set(GRB_IntParam_OutputFlag, 0);
		model.optimize();

		// Get information
		*xObjValue= model.get(GRB_DoubleAttr_ObjVal);
		endTime = clock();
		//if (xIns.mState > 2) { xIns.mCpuTime = STP_SOVLERTIME; }
		//else { xIns.mCpuTime = double(double(endTime - startTime) / CLOCKS_PER_SEC); }
	}
	catch (GRBException ee) {
		cout << ee.getErrorCode() << endl;
		cout << ee.getMessage() << endl;

		xIns.mState = ee.getErrorCode();
		//endTime = clock();
		//xIns.mCpuTime = double(double(endTime - startTime) / CLOCKS_PER_SEC);
		cout << "This is Model P" << endl;
		system("pause");
	}
	catch (...) {
		cout << "Some problems happen" << endl;
		cout << "This is Model P" << endl;
		system("pause");
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



