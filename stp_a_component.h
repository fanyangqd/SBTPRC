#pragma once
#ifndef STP_A_COMPONENT_H
#define STP_A_COMPONENT_H

// Required include's 
#include "gurobi_c++.h"
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <queue>
#include <vector>
#include <iterator>
#include <numeric>
#include <ctime>
#include <bitset>
#include <unordered_map>
#include <set>
using namespace std;

#define STP_SOVLERTIME 5*STP_N
#define STP_N 60          //Component number;
#define	STP_RESCAP 10	  //Resouce capacity;
#define STP_RES_LB 2	  //Resource requirement range
#define STP_RES_UB 4
#define STP_QRANGE_LB 0.01 //Range for q
#define STP_QRANGE_UB 0.5
#define STP_FIXEDCOST 8*STP_N
#define STP_VMAX 2// upper bound for visiting each level of batch b


#define STP_EPSILON 0.00001





struct cmp {
	bool operator()(vector<double>& a, vector<double>& b);
};
typedef priority_queue<vector<double>, vector<vector<double>>, cmp> allBatches;
//Definition for components
class StpComponent {
public:
	int	mRes;		//resoruce requirement
	double mCost;		//cost
	int mWeight;	//weight
	double mProb;	//probability
	double mFixedCost;  // Fixed cost
	
	//functions:
	StpComponent();
};


//Definition for instances
class StpInstance {
public:
	StpComponent mComp[STP_N+1];
	double mObjVal;		    // Objective value
	double mCpuTime;		// Computation time
	double mGap;		    // Gap between upper-bound and lower-bound
	int mState;				// solution state: Optimal/feasible/infeasible
	double tempTTB;         // temporary time to best
	double tempTTB02;
	double tempTTB03;
	double timeToBest;   // time to Best
	clock_t recordStart;
	allBatches insAllBa; //record all the batches;

	//functions:
	StpInstance();
	void Initialization();  // Generate job data randomly
	void Reset();			// Set all job data as zero
	void Output(int Instance_Number); // Output job data of an instance to a txt file 
	void Input(string Instance_Name); // Input job data from a txt file to an instance
};

class StpBatch {
public:
	int batchID;
	int compNumber;
	int batchCost;
	double batchProb;
	double batchRatio;
	StpBatch();
	vector<int> comps;
};

double max_number(double a, double b);
double min_number(double a, double b);
int min_number(int a, int b);

void compute_ObjVal(double& xObj, vector<vector<double>>& xBatchSeq);
bool lessSort(vector<double>& a, vector<double>& b);

#endif // 
