#pragma once
#ifndef STP_TABU_SEARCH_H
#define STP_TABU_SEARCH_H

#include "stp_a_component.h"


//bool LessSort(vector<double>& a, vector<double>& b);
int weighed_random(vector<int>& xVec);


class TabuSolver {
public:
  
    vector<vector<double>> bestSol;// best Solution
    vector<vector<double>> currSol;// current Solution
    double bestObj;
    double currObj;
    vector<vector<double>> final_bestSol;
    double final_bestObj;

    int batchLB; //minmum number for batching all of components
    int tabuList[STP_N + 1];
    int compMoved; //record which component is moved
    int repairMoved;
    int tabuTenure;
    int currIteration;
    int tabuStop; // stopping criterion 
    int counterB;
    int visit[STP_N+1]; //record visiting history for each level of batch b;
    StpInstance *ins;

    clock_t startTime, endTime;

    // For hybrid proximity
    int mMode; // 0: Original TS; 1: initial sol for PB; 2: medial in PB

    GRBVar** mVar;
    GRBLinExpr* mDelta;
    double* mIncumbent; // Obtained objective value from proximity search
    double** mSol;
    bool* mImproved;
    

    //for re-initialization
    int reIter;
   // vector<vector<vector<double>>> reInitialSol;

    //Functions:
    void initial_solution(StpInstance* xIns);
   // void compute_ObjVal(double &xObj, vector<vector<double>> &xBatchSeq);
    TabuSolver(StpInstance *xIns);
    TabuSolver(StpInstance* xIns,  GRBVar** xVar,  GRBLinExpr* xObjFunc, 
        double* xIncumbent, int xMode, double** xSol, bool* xImproved, clock_t xStart);
   

    void solve_re_initial();
    void get_best_neighbor();

    void comp_insertion(int xComp, int xOriBatch,
        int xDestBatch, vector<vector<double>>& xBatchSeq);

    void basic_tabu_search();
    bool batch_reduction();

};

int createArcs(StpInstance* xIns, vector<vector<int> >& arcs, vector<bool>& souvenir);
int transformSolution(vector<vector<vector<int> > >& arcsUsed, vector<vector<int> >& sol);
int arcflow(StpInstance* xIns, vector<vector<int> >& sol); // arc flow model for bin packing
void binPacking(StpInstance* xIns); // classical model for bin packing

// Hybrid Proximity
void hybrid_proximity(StpInstance& xIns);
void refine(StpInstance& xIns, GRBVar** xVar, double* xObjValue, clock_t xStart);

#endif //#pragma once
