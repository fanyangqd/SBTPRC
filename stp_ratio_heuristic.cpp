#include"stp_ratio_heuristic.h"
#include"stp_tabu_search.h"


void generate_batches(StpInstance *xIns, int xCounter, vector<bool>& xMySet,  
	allBatches& xAB, int xTotalRes)
{
	if (xCounter > STP_N) {
		// shows the subset
		vector<double> tempBatch = { 0.0,0.0,1.0,STP_FIXEDCOST };
		for (int i = 1; i <=STP_N; i++) {
			if (xMySet[i] == true) {
				//cout << i << " ";
				tempBatch.push_back(i);
				tempBatch[0] += xIns->mComp[i].mRes;
				tempBatch[2] = tempBatch[2] * xIns->mComp[i].mProb;
				tempBatch[3] += xIns->mComp[i].mCost;
			}
		}
		//cout << endl;
		if (tempBatch.size() > 4) { 
			tempBatch[1] = tempBatch[3] / (1 - tempBatch[2]);
			xAB.emplace(tempBatch); }
	}
	else {
		if (xTotalRes + xIns->mComp[xCounter].mRes> STP_RESCAP) {
			xMySet[xCounter] = false;
			generate_batches(xIns, xCounter + 1, xMySet, xAB, xTotalRes);
		}
		else {
			xMySet[xCounter] = true;
			xTotalRes += xIns->mComp[xCounter].mRes;
			generate_batches(xIns, xCounter + 1, xMySet, xAB, xTotalRes);
			xTotalRes -= xIns->mComp[xCounter].mRes;
			xMySet[xCounter] = false; // backtracking
			generate_batches(xIns, xCounter + 1, xMySet, xAB, xTotalRes);
		}
	}

	return;
}


void ratio_heuristic(StpInstance* xIns){
	clock_t startTime, endTime;
	startTime = clock();
    
	//generate all of batches
	vector<bool> mySet; mySet.resize(STP_N + 1);
	allBatches allBa;
	generate_batches(xIns, 1, mySet, allBa, 0);
	xIns->insAllBa = allBa;

	vector<vector<double>> bestSol, curSol;
	double bestObj = INT_MAX;
	//double curObj = 0.0;
	int record[STP_N + 1];
	memset(record, 0, sizeof(int) * (STP_N + 1));
	while (!allBa.empty()) {
		vector<double> tempBatch = allBa.top();
		curSol.push_back(tempBatch);
		allBa.pop();
		for (int i = 4; i < tempBatch.size(); i++){
			record[int(tempBatch[i])] = 1;
		}

		// add singleton test
		vector<vector<double>> copySol = curSol;
		double copyObj = 0.0;
		for (int i = 1; i <=STP_N; i++){
			if (record[i] < 1) {
				vector<double> singletonBatch = { double(xIns->mComp[i].mRes),
					0.0,xIns->mComp[i].mProb,double(STP_FIXEDCOST+xIns->mComp[i].mCost),double(i) };
				singletonBatch[1] = singletonBatch[3] / (1 - singletonBatch[2]);
				copySol.push_back(singletonBatch);
			}
		}
		compute_ObjVal(copyObj, copySol);
		//if (copyObj < bestObj) {
		if(bestObj-copyObj>STP_EPSILON){
			bestObj = copyObj;
			bestSol = copySol;
			xIns->tempTTB = double(double(clock() - startTime) / CLOCKS_PER_SEC);
		}

		//delete batches in set, which intersects with batches in solution
		tempBatch.clear();
		allBatches tempAllBa;
		while (!allBa.empty()) {
			tempBatch = allBa.top();
			allBa.pop();
			bool intersected = false;

			for (int i = 4; i < tempBatch.size(); i++){
				for (int t = 0; t < curSol.size(); t++){
					vector<double>::iterator it = find(curSol[t].begin()+4, curSol[t].end(), tempBatch[i]);
					if (it != curSol[t].end()) {
						intersected = true;
						break;
					}
				}
				if (intersected) { break; }
			}
			if (!intersected) {
				tempAllBa.emplace(tempBatch);
			}
		}

		allBa = tempAllBa;
	}
	cout << "Ratio heuristic: " << bestObj << endl;
	xIns->mObjVal = bestObj;
	endTime = clock();
	xIns->mCpuTime=double(double(endTime - startTime) / CLOCKS_PER_SEC); 
	xIns->timeToBest = xIns->tempTTB;
	return;
}


double ratio_heuristic(StpInstance* xIns, vector<vector<double>> &xSol) {
	clock_t startTime;// , endTime;
	startTime = clock();

	//generate all of batches
	vector<bool> mySet; mySet.resize(STP_N + 1);
	allBatches allBa;
	generate_batches(xIns, 1, mySet, allBa, 0);

	vector<vector<double>> bestSol, curSol;
	double bestObj = INT_MAX;
	//double curObj = 0.0;
	int record[STP_N + 1];
	memset(record, 0, sizeof(int) * (STP_N + 1));
	while (!allBa.empty()) {
		vector<double> tempBatch = allBa.top();
		curSol.push_back(tempBatch);
		allBa.pop();
		for (int i = 4; i < tempBatch.size(); i++) {
			record[int(tempBatch[i])] = 1;
		}

		// add singleton test
		vector<vector<double>> copySol = curSol;
		double copyObj = 0.0;
		for (int i = 1; i <= STP_N; i++) {
			if (record[i] < 1) {
				vector<double> singletonBatch = { double(xIns->mComp[i].mRes),
					0.0,xIns->mComp[i].mProb,double(STP_FIXEDCOST + xIns->mComp[i].mCost),double(i) };
				singletonBatch[1] = singletonBatch[3] / (1 - singletonBatch[2]);
				
				copySol.push_back(singletonBatch);
			}
		}
		

		compute_ObjVal(copyObj, copySol);
		
	//	if (copyObj < bestObj) {
		if (bestObj - copyObj > STP_EPSILON) {
			bestObj = copyObj;
			bestSol = copySol;
			//xIns->tempTTB = double(double(clock() - xIns->recordStart) / CLOCKS_PER_SEC);
		}

		//delete batches in set, which intersects with batches in solution
		tempBatch.clear();
		allBatches tempAllBa;
		while (!allBa.empty()) {
			tempBatch = allBa.top();
			allBa.pop();
			bool intersected = false;

			for (int i = 4; i < tempBatch.size(); i++) {
				for (int t = 0; t < curSol.size(); t++) {
					vector<double>::iterator it = find(curSol[t].begin() + 4, curSol[t].end(), tempBatch[i]);
					if (it != curSol[t].end()) {
						intersected = true;
						break;
					}
				}
				if (intersected) { break; }
			}
			if (!intersected) {
				tempAllBa.emplace(tempBatch);
			}
		}

		allBa = tempAllBa;
	}
	xSol = bestSol;
	
	//xIns->mObjVal = bestObj;
	//endTime = clock();
	//xIns->mCpuTime = double(double(endTime - startTime) / CLOCKS_PER_SEC);
	return bestObj;
}


// For TS
void random_initialization(StpInstance* xIns, vector<vector<double>>& xSol) {
	clock_t startTime;// , endTime;
	startTime = clock();

	//generate all of batches
	vector<bool> mySet; mySet.resize(STP_N + 1);
	allBatches allBa;
	generate_batches(xIns, 1, mySet, allBa, 0);

	vector<vector<double>> bestSol, curSol;
	double bestObj = INT_MAX;

	int record[STP_N + 1];
	memset(record, 0, sizeof(int) * (STP_N + 1));
	int firstBatch = 1;
	while (!allBa.empty()) {
		vector<double> tempBatch;
		if (firstBatch) {
			firstBatch = 0;
			int totalSize = int(allBa.size());// return the size for all of batches
			std::random_device rd;
			std::mt19937 gen(rd());
			std::uniform_int_distribution<> dis(1, totalSize);
			int randomNumber = dis(gen); //generate the batch number randomly
			int counter = 0;
			allBatches tempContainer;

			while (counter < randomNumber) {
				counter++;
				if (counter < randomNumber) { tempContainer.push(allBa.top()); }
				else { tempBatch = allBa.top(); }
				allBa.pop();
			}
			while (!tempContainer.empty()) {
				allBa.push(tempContainer.top());
				tempContainer.pop();
			}
		}
		else {
			tempBatch = allBa.top();
			allBa.pop();
		}

		curSol.push_back(tempBatch);

		for (int i = 4; i < tempBatch.size(); i++) {
			record[int(tempBatch[i])] = 1;
		}

		// add singleton test
		vector<vector<double>> copySol = curSol;
		double copyObj = 0.0;
		for (int i = 1; i <= STP_N; i++) {
			if (record[i] < 1) {
				vector<double> singletonBatch = { double(xIns->mComp[i].mRes),
					0.0,xIns->mComp[i].mProb,double(STP_FIXEDCOST + xIns->mComp[i].mCost),double(i) };
				singletonBatch[1] = singletonBatch[3] / (1 - singletonBatch[2]);
				copySol.push_back(singletonBatch);
			}
		}
		compute_ObjVal(copyObj, copySol);
	//	if (copyObj < bestObj) {
		if (bestObj - copyObj > STP_EPSILON) {
			bestObj = copyObj;
			bestSol = copySol;	
		}

	
		//delete batches in set, which intersects with batches in solution
		tempBatch.clear();
		allBatches tempAllBa;
		while (!allBa.empty()) {
			tempBatch = allBa.top();
			allBa.pop();
			bool intersected = false;

			for (int i = 4; i < tempBatch.size(); i++) {
				for (int t = 0; t < curSol.size(); t++) {
					vector<double>::iterator it = find(curSol[t].begin() + 4, curSol[t].end(), tempBatch[i]);
					if (it != curSol[t].end()) {
						intersected = true;
						break;
					}
				}
				if (intersected) { break; }
			}
			if (!intersected) {
				tempAllBa.emplace(tempBatch);
			}
		}

		allBa = tempAllBa;
	}
	xSol = bestSol;
	return;
}

void adapted_ratio_heuristic(StpInstance* xIns) {
	clock_t startTime, endTime;
	startTime = clock();
	xIns->recordStart = startTime;
	//generate all of batches
	vector<bool> mySet; mySet.resize(STP_N + 1);
	allBatches allBatemp;
	generate_batches(xIns, 1, mySet, allBatemp, 0);

	vector<vector<double>> bestSol;
	double bestObj = ratio_heuristic(xIns, bestSol);
	xIns->tempTTB = double(double(clock() - startTime) / CLOCKS_PER_SEC);
	//int round = 0;
	//vector<double> objVector;
	//objVector.push_back(bestObj);
	while (double(double(clock() - startTime) / CLOCKS_PER_SEC) < STP_SOVLERTIME) {
		//round++;
		allBatches allBa = allBatemp;
		vector<vector<double>> curSol;
		int record[STP_N + 1];
		memset(record, 0, sizeof(int) * (STP_N + 1));
		//bestObj = INT_MAX;
		int firstBatch = 1;
		while (!allBa.empty()) {

			vector<double> tempBatch;
			if (firstBatch) {
				firstBatch = 0;
				int totalSize = int(allBa.size());// return the size for all of batches
				std::random_device rd;
				std::mt19937 gen(rd());
				std::uniform_int_distribution<> dis(1, totalSize);
				int randomNumber = dis(gen); //generate the batch number randomly
				int counter = 0;
				allBatches tempContainer;
			
				while (counter < randomNumber) {
					counter++;
					if (counter < randomNumber) { tempContainer.push(allBa.top()); }
					else { tempBatch = allBa.top(); }
					allBa.pop();
				}
				while (!tempContainer.empty()) {
					allBa.push(tempContainer.top());
					tempContainer.pop();
				}
			}
			else {
				tempBatch = allBa.top();
				allBa.pop();
			}
		
			curSol.push_back(tempBatch);

			for (int i = 4; i < tempBatch.size(); i++) {
				record[int(tempBatch[i])] = 1;
			}

			// add singleton test
			vector<vector<double>> copySol = curSol;
			double copyObj = 0.0;
			for (int i = 1; i <= STP_N; i++) {
				if (record[i] < 1) {
					vector<double> singletonBatch = { double(xIns->mComp[i].mRes),
						0.0,xIns->mComp[i].mProb,double(STP_FIXEDCOST + xIns->mComp[i].mCost),double(i) };
					singletonBatch[1] = singletonBatch[3] / (1 - singletonBatch[2]);
					copySol.push_back(singletonBatch);
				}
			}
			compute_ObjVal(copyObj, copySol);
		//	if (copyObj < bestObj) {
			if (bestObj - copyObj > STP_EPSILON) {
				bestObj = copyObj;
				bestSol = copySol;
				xIns->tempTTB = double(double(clock() - startTime) / CLOCKS_PER_SEC);
			}

			//delete batches in set, which intersects with batches in solution
			tempBatch.clear();
			allBatches tempAllBa;
			while (!allBa.empty()) {
				tempBatch = allBa.top();
				allBa.pop();
				bool intersected = false;

				for (int i = 4; i < tempBatch.size(); i++) {
					for (int t = 0; t < curSol.size(); t++) {
						vector<double>::iterator it = find(curSol[t].begin() + 4, curSol[t].end(), tempBatch[i]);
						if (it != curSol[t].end()) {
							intersected = true;
							break;
						}
					}
					if (intersected) { break; }
				}
				if (!intersected) {
					tempAllBa.emplace(tempBatch);
				}
			}

			allBa = tempAllBa;
		}


		//vector<double>::iterator it = find(objVector.begin(), objVector.end(), bestObj);
		//if (it == objVector.end()) {objVector.push_back(bestObj);}

	}
	//cout << "Ratio heuristic: " << bestObj << endl;

	
	//cout << "Round: " << round << " Obj number: " << objVector.size() << endl;
	//xIns->mObjVal = round;
	//xIns->mCpuTime = objVector.size();
	xIns->mObjVal = bestObj;
	endTime = clock();
	xIns->mCpuTime = double(double(endTime - startTime) / CLOCKS_PER_SEC);
	xIns->timeToBest = xIns->tempTTB;
	return;
}
