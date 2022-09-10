#include "checker.h"
#include "pb/util/unumeric.h"

Checker::Checker(const DetProb &prob, double epsilon): prob(prob), epsilon(epsilon){
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
    if (isLess(Ax(i), prob.bl(i), epsilon)) return LbConstraint;
    if (isGreater(Ax(i), prob.bu(i), epsilon)) return UbConstraint;
  }
  for (int i = 0; i < n; i ++){
    if (isLess(sol(i), prob.l(i), epsilon)) return LbVariable;
    if (isGreater(sol(i), prob.u(i), epsilon)) return UbVariable;
  }
  return Feasibility;
}

int Checker::checkIlpFeasibility(const VectorXd &sol){
  if (sol.size() != n) return Unsolved;
  for (int i = 0; i < n; i ++){
    if (!isInteger(sol(i), epsilon)){
      return Integrality;
    }
  }
  return checkLpFeasibility(sol);
}