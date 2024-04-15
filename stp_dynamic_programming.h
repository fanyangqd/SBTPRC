#pragma once
#ifndef STP_DYNAMIC_PROGRAMMING_H
#define STP_DYNAMIC_PROGRAMMING_H

#include "stp_a_component.h"


//Define bitsets
typedef bitset<STP_N> bits;

//first: set of past activities, second: value function
typedef unordered_map<bits, double, hash<bits>> stageset;


//void state_initialize(StpInstance& xIns, int xCompNum, vector<int>& xBatchSet,
	//vector<stageset*>& xAllStates);// , vector<int>& xNrInStage);

void state_initialize(StpInstance& xIns, int xCompNum, int xBatchSet[],
	vector<stageset*>& xAllStates);

//void batch_initialize(StpInstance& xIns, int xCompNum, vector<int>& xBatchSet, int xTotalRes,
	//unordered_map<bits, double, hash<bits>>& xBatchProb);

//void state_value_compute(bits& xState, double& xStateValue, vector<stageset*>& xAllStates, StpInstance& xIns, int xCompNum,
	//vector<int>& xBatchSet, int xTotalRes, unordered_map<bits, double, hash<bits>>& xBatchProb);

void dynamic_programming(StpInstance& xIns);

void batch_initialize(StpInstance& xIns, int xCompNum, int xBatchSet[], int xTotalRes,
	unordered_map<bits, double, hash<bits>>& xBatchProb, unordered_map<bits, double, hash<bits>>& xBatchCost);

void state_value_compute(bits& xState, double& xStateValue, vector<stageset*>& xAllStates, StpInstance& xIns, int xCompNum,
	int xBatchSet[], int xTotalRes, unordered_map<bits, double, hash<bits>>& xBatchProb, unordered_map<bits, double, hash<bits>>& xBatchCost);

#endif //
