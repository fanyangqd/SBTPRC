#include "stp_dynamic_programming.h"

void state_initialize(StpInstance& xIns,  int xCompNum, int xBatchSet[], 
	vector<stageset*> &xAllStates){//, vector<int> &xNrInStage) {
	//Generate all the states for each stage
	if (xCompNum > STP_N) {
		// shows the subset
		//vector<int> cc;
		bits tempState;
		for (int i = 1; i <= STP_N; i++) {
			if (xBatchSet[i] == 1) {
				//cc.push_back(i);
				tempState[i - 1] = 1;
			}
		}

		int tempSize = int(tempState.count());
		if (tempSize > 0) {
			//Search for state predP in prev_stage
			stageset::iterator existing_state = xAllStates[STP_N-tempSize+1]->find(tempState);

			//If it doesn't exist already, add it.
			if (existing_state == xAllStates[STP_N - tempSize + 1]->end()) {
				xAllStates[STP_N - tempSize + 1]->insert({ tempState, 0 });
				//xNrInStage[STP_N - tempSize + 1]++;
			}//end if
		}

	}
	else{
			xBatchSet[xCompNum] = 1;
			state_initialize(xIns, xCompNum + 1, xBatchSet, xAllStates);// , xNrInStage);
			xBatchSet[xCompNum] = 0; // backtracking
			state_initialize(xIns, xCompNum + 1, xBatchSet, xAllStates);// , xNrInStage);
	}

	return;
}
/*
void batch_initialize(StpInstance& xIns, int xCompNum, vector<int>& xBatchSet, int xTotalRes,
	unordered_map<bits, double, hash<bits>> &xBatchProb) {
	//Generate all the batch
	if (xCompNum > STP_N) {
		// shows the subset
		double cost = xIns.mComp[1].mFixedCost;
		double prob = 1;
		bits tempState;
		for (int i = 1; i <= STP_N; i++) {
			if (xBatchSet[i] == 1) {
				tempState[i - 1] = 1;
				cost += xIns.mComp[i].mCost;
				prob = prob * xIns.mComp[i].mProb;
			}
		}

		if (tempState.count() > 0) {
			//Search for state predP in prev_stage
			unordered_map<bits, double, hash<bits>>::iterator existing_state = xBatchProb.find(tempState);

			//If it doesn't exist already, add it.
			if (existing_state == xBatchProb.end()) {
				xBatchProb.insert({ tempState, cost+prob });
				//xBatchProb.insert({ tempState, prob });
			}//end if
		}

	}
	else {
		if (xTotalRes + xIns.mComp[xCompNum].mRes > STP_RESCAP) {
			xBatchSet[xCompNum] = 0;
			batch_initialize(xIns, xCompNum+1, xBatchSet, xTotalRes, xBatchProb);
		}
		else {
			xBatchSet[xCompNum] = 1;
			xTotalRes += xIns.mComp[xCompNum].mRes;
			batch_initialize(xIns, xCompNum+1, xBatchSet, xTotalRes,  xBatchProb);
			xTotalRes -= xIns.mComp[xCompNum].mRes;
			xBatchSet[xCompNum] = 0; // backtracking
			batch_initialize(xIns, xCompNum+1, xBatchSet, xTotalRes,  xBatchProb);
		}
	}

	return;
}


void state_value_compute(bits& xState, double& xStateValue, vector<stageset*>& xAllStates,StpInstance& xIns, int xCompNum,
	vector<int>& xBatchSet, int xTotalRes, unordered_map<bits, double, hash<bits>>& xBatchProb) {
	
	if (xCompNum > STP_N) {
		bits batchState;
		for (int i = 1; i <= STP_N; i++) {
			if (xBatchSet[i] == 1) {
				batchState[i - 1] = 1;
			}
		}

		if (batchState.count() > 0) {
			bits remainState = (~batchState) & xState;
			double batchCost, batchProb;
			double remainValue;
			unordered_map<bits, double, hash<bits>>::iterator existing_batch = xBatchProb.find(batchState);

			if (existing_batch != xBatchProb.end()) {
				batchCost = int(existing_batch->second);
				batchProb = existing_batch->second - batchCost;
			}


			stageset::iterator remain_state = xAllStates[STP_N - int(remainState.count()) + 1]->find(remainState);
			if (remain_state != xAllStates[STP_N - int(remainState.count()) + 1]->end()) {
				remainValue = remain_state->second;
			}

			if (batchCost + batchProb * remainValue < xStateValue) {
				xStateValue = batchCost + batchProb * remainValue;
			}

		}

	}
	else {
		//cout << xState[xCompNum - 1] << endl;

		if (xTotalRes + xIns.mComp[xCompNum].mRes > STP_RESCAP || xState[xCompNum - 1] == 0) {

			xBatchSet[xCompNum] = 0;
			state_value_compute(xState,xStateValue,xAllStates,xIns,xCompNum+1,xBatchSet,xTotalRes,xBatchProb);
		}
		else {
			
			xBatchSet[xCompNum] = 1;
			xTotalRes += xIns.mComp[xCompNum].mRes;
			
			state_value_compute(xState, xStateValue, xAllStates, xIns, xCompNum + 1, xBatchSet, xTotalRes, xBatchProb);

			xTotalRes -= xIns.mComp[xCompNum].mRes;
			xBatchSet[xCompNum] = 0; // backtracking
			
			state_value_compute(xState, xStateValue, xAllStates, xIns, xCompNum + 1, xBatchSet, xTotalRes, xBatchProb);
		}
	}
	return;
}
*/
void dynamic_programming(StpInstance& xIns) {
	clock_t startTime, endTime;
	startTime = clock();
	//Create stage vector
	vector <stageset*> stage(STP_N+2);	 //vector with stage[k] pointing to hash table, from 1 to STP_N

	//vector <int> nr_of_states_in_stage(STP_N+1);	//vector with element k containing the number of states in stage k

	//initialize stage vector and nr states vector
	for (int k = STP_N+1; k > 0; --k) {
		stage[k] = new stageset;
		//nr_of_states_in_stage[k] = 0;
	}
	// boundary condition
	bits boundaryCondition;
	stage[STP_N + 1]->insert({ boundaryCondition, 0 });

	int my_set[STP_N + 1] = { 0 };
	state_initialize(xIns, 1, my_set, stage);//, nr_of_states_in_stage);
	

	unordered_map<bits, double, hash<bits>> batchProb;
	unordered_map<bits, double, hash<bits>> batchCost;
	memset(my_set, 0, sizeof(my_set));
	
	batch_initialize(xIns, 1, my_set, 0, batchProb,batchCost);
	

	// Compute the value of each state
	for (int k = STP_N; k > 0; --k) {
		if (double(double(clock() - startTime) / CLOCKS_PER_SEC) - STP_SOVLERTIME > STP_EPSILON) { break; }
		for (auto stateP = stage[k]->begin(); stateP != stage[k]->end(); ++stateP){
			if (double(double(clock() - startTime) / CLOCKS_PER_SEC) - STP_SOVLERTIME > STP_EPSILON) { break; }
			bits currentState = stateP->first;
			//my_set.resize(STP_N + 1);
			memset(my_set, 0, sizeof(my_set));
			double currentValue = INT_MAX;
			state_value_compute(currentState, currentValue, stage,xIns, 1, my_set, 0, batchProb,batchCost );

			stateP->second = currentValue;
			//cout << currentState << "  " << currentValue <<  endl;
		
		}
	}
	endTime = clock();
	xIns.mCpuTime = double(double(endTime - startTime) / CLOCKS_PER_SEC);
	bits finalState; finalState.set();
	stageset::iterator ff = stage[1]->find(finalState);
	xIns.mObjVal = ff->second;
	for (int k = STP_N+1; k > 0; --k) {delete stage[k];}
	return;
}



void batch_initialize(StpInstance& xIns, int xCompNum, int xBatchSet[], int xTotalRes,
	unordered_map<bits, double, hash<bits>>& xBatchProb, unordered_map<bits, double, hash<bits>>& xBatchCost) {
	//Generate all the batch
	if (xCompNum > STP_N) {
		// shows the subset
		double cost = xIns.mComp[1].mFixedCost;
		double prob = 1.0;
		bits tempState;
		for (int i = 1; i <= STP_N; i++) {
			if (xBatchSet[i] == 1) {
				tempState[i - 1] = 1;
				cost += xIns.mComp[i].mCost;
				prob = prob * xIns.mComp[i].mProb;
			}
		}

		if (tempState.count() > 0) {
			//Search for state predP in prev_stage
			unordered_map<bits, double, hash<bits>>::iterator existing_state = xBatchProb.find(tempState);

			//If it doesn't exist already, add it.
			if (existing_state == xBatchProb.end()) {
				xBatchProb.insert({ tempState,  prob });
			}//end if

			unordered_map<bits, double, hash<bits>>::iterator existing_state02 = xBatchCost.find(tempState);
			//If it doesn't exist already, add it.
			if (existing_state02 == xBatchCost.end()) {
				xBatchCost.insert({ tempState,  cost });
			}//end if
		}

	}
	else {
		if (xTotalRes + xIns.mComp[xCompNum].mRes > STP_RESCAP) {
			xBatchSet[xCompNum] = 0;
			batch_initialize(xIns, xCompNum + 1, xBatchSet, xTotalRes, xBatchProb, xBatchCost);
		}
		else {
			xBatchSet[xCompNum] = 1;
			xTotalRes += xIns.mComp[xCompNum].mRes;
			batch_initialize(xIns, xCompNum + 1, xBatchSet, xTotalRes, xBatchProb, xBatchCost);
			xTotalRes -= xIns.mComp[xCompNum].mRes;
			xBatchSet[xCompNum] = 0; // backtracking
			batch_initialize(xIns, xCompNum + 1, xBatchSet, xTotalRes, xBatchProb, xBatchCost);
		}
	}

	return;
}


void state_value_compute(bits& xState, double& xStateValue, vector<stageset*>& xAllStates, StpInstance& xIns, int xCompNum,
	int xBatchSet[], int xTotalRes, unordered_map<bits, double, hash<bits>>& xBatchProb, 
	unordered_map<bits, double, hash<bits>>& xBatchCost) {

	if (xCompNum > STP_N) {
		bits batchState;
		for (int i = 1; i <= STP_N; i++) {
			if (xBatchSet[i] == 1) {
				batchState[i - 1] = 1;
			}
		}

		if (batchState.count() > 0) {
			bits remainState = (~batchState) & xState;
			double batchCost, batchProb;
			double remainValue;
			unordered_map<bits, double, hash<bits>>::iterator existing_batch = xBatchProb.find(batchState);
			if (existing_batch != xBatchProb.end()) {
				batchProb = existing_batch->second;
			}

			unordered_map<bits, double, hash<bits>>::iterator existing_batch02 = xBatchCost.find(batchState);
			if (existing_batch02 != xBatchCost.end()) {
				batchCost = existing_batch02->second;
			}


			stageset::iterator remain_state = xAllStates[STP_N - int(remainState.count()) + 1]->find(remainState);
			if (remain_state != xAllStates[STP_N - int(remainState.count()) + 1]->end()) {
				remainValue = remain_state->second;
			}

			if (batchCost + batchProb * remainValue < xStateValue) {
				xStateValue = batchCost + batchProb * remainValue;
			}

		}

	}
	else {
		//cout << xState[xCompNum - 1] << endl;

		if (xTotalRes + xIns.mComp[xCompNum].mRes > STP_RESCAP || xState[xCompNum - 1] == 0) {

			xBatchSet[xCompNum] = 0;
			state_value_compute(xState, xStateValue, xAllStates, xIns, xCompNum + 1, xBatchSet, xTotalRes, xBatchProb,xBatchCost);
		}
		else {

			xBatchSet[xCompNum] = 1;
			xTotalRes += xIns.mComp[xCompNum].mRes;

			state_value_compute(xState, xStateValue, xAllStates, xIns, xCompNum + 1, xBatchSet, xTotalRes, xBatchProb,xBatchCost);

			xTotalRes -= xIns.mComp[xCompNum].mRes;
			xBatchSet[xCompNum] = 0; // backtracking

			state_value_compute(xState, xStateValue, xAllStates, xIns, xCompNum + 1, xBatchSet, xTotalRes, xBatchProb, xBatchCost);
		}
	}

	return;
}


