#include "checker.h"
#include "pb/util/unumeric.h"

double Checker::kEpsilon = 1e-10;

Checker::Checker(const DetProb &prob): prob(prob){
  n = (int) prob.c.size();
  m = (int) prob.bl.size();
}

double Checker::getScore(const VectorXd &sol){
  return prob.c.dot(sol);
}

int Checker::checkLpFeasibility(const VectorXd &sol){
  if (sol.size() != n) return Unsolved;
  VectorXd Ax = prob.A * sol;
  for (int i = 0; i < m; i ++){
    if (isLess(Ax(i), prob.bl(i), kEpsilon)) return LbConstraint;
    if (isGreater(Ax(i), prob.bu(i), kEpsilon)) return UbConstraint;
  }
  for (int i = 0; i < n; i ++){
    if (isLess(sol(i), prob.l(i), kEpsilon)) return LbVariable;
    if (isGreater(sol(i), prob.u(i), kEpsilon)) return UbVariable;
  }
  return Feasibility;
}

int Checker::checkIlpFeasibility(const VectorXd &sol){
  if (sol.size() != n) return Unsolved;
  for (int i = 0; i < n; i ++){
    if (!isInteger(sol(i), kEpsilon)) return Integrality;
  }
  return checkLpFeasibility(sol);
}