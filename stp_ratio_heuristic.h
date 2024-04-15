#pragma once
#ifndef STP_RATIO_HEURISTIC_H
#define STP_RATIO_HEURISTIC_H

#include "stp_a_component.h"



//void compute_Obj(double& xObj, vector<vector<double>>& xBatchSeq);
void generate_batches(StpInstance* xIns, int xCounter, vector<bool>& xMySet,
	allBatches& xAB, int xTotalRes);

void generate_batches(StpInstance* xIns, int xCounter, vector<bool>& xMySet,
	allBatches & xAB, int xTotalRes);

void ratio_heuristic(StpInstance* xIns);
// For  initialization of assignmend-based model
double ratio_heuristic(StpInstance* xIns, vector<vector<double>>& xSol);

// For tabu search re-initialization
void random_initialization(StpInstance* xIns, vector<vector<double>>& xSol);

void adapted_ratio_heuristic(StpInstance* xIns);

#endif //

