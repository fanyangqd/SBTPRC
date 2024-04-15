#pragma once
#ifndef STP_MODEL_SOLVER_H
#define STP_MODEL_SOLVER_H

#include "stp_a_component.h"

void compute_obj_val(StpInstance& xIns,GRBVar **xVar);
void compute_obj_val(StpInstance& xIns, GRBVar** xVar, double& xObj);
void assignment_based_model(StpInstance &xIns);
//void assignment_based_model_validCut(StpInstance& xIns);


#endif //

