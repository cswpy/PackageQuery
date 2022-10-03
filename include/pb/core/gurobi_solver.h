#pragma once

#include "pb/det/det_prob.h"
#include "pb/util/udeclare.h"
#include "gurobi_c++.h"
#include "checker.h"

using namespace pb;

class GurobiSolver: public Checker{
private:
  double exe_init, exe_lp, exe_ilp;
  GRBenv *env;
  GRBmodel *model;
public:
  int lp_status, ilp_status;
  VectorXd lp_sol, ilp_sol;
  double lp_score, ilp_score;
  int iteration_count;
public:
  ~GurobiSolver();
  GurobiSolver(const DetProb &prob, bool with_objective=true);
  void solveIlp(double mipGap=1e-4, double time_limit=-1.0);
  void solveLp();
  double getLpTime();
  double getIlpTime();
  double getScore(const VectorXd &sol);
  int checkLpFeasibility(const VectorXd &sol);
  int checkIlpFeasibility(const VectorXd &sol);
  void writeModel(string file_name);
  bool hasIlpSolution();
};