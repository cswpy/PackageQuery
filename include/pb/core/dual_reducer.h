#pragma once

#include "pb/util/udeclare.h"
#include "pb/core/dual.h"
#include "pb/det/det_prob.h"

using namespace pb;

class DualReducer{
public:
  static int kIlpSize;
  static double kEpsilon, kTimeLimit;
  VectorXd ilp_sol, lp_sol;
  double ilp_score, lp_score, exe_ilp, exe_lp;
  int status;
public:
  ~DualReducer();
  DualReducer(int core, const DetProb &prob);
  DualReducer(int core, const DetProb &prob, VectorXd oracle);
};
