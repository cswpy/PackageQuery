#include "checker.h"
#include "pb/util/unumeric.h"

Checker::Checker(const DetProb &prob, double epsilon): prob(prob), epsilon(epsilon){
  n = (int) prob.c.size();
  m = (int) prob.bl.size();
  if (n == (int) prob.ids.size()){
    id_map.reserve(n);
    for (int i = 0; i < n; i ++) id_map[prob.ids[i]] = i;
  }
}

double Checker::getScore(const VectorXd &sol){
  return prob.c.dot(sol);
}

int Checker::checkLpFeasibility(const VectorXd &sol){
  if (sol.size() != n) return Unsolved;
  VectorXd Ax = prob.A * sol;
  for (int i = 0; i < m; i ++){
    if (prob.bl(i) != -DBL_MAX && isLess(Ax(i), prob.bl(i), epsilon)) return LbConstraint;
    if (prob.bu(i) != DBL_MAX && isGreater(Ax(i), prob.bu(i), epsilon)) return UbConstraint;
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

double Checker::getScore(map<long long, double> &sol){
  VectorXd vsol (n); vsol.fill(0);
  for (const auto& [i, v] : sol) vsol(id_map[i]) = v;
  return getScore(vsol);
}

double Checker::getScore(map<long long, long long> &sol){
  VectorXd vsol (n); vsol.fill(0);
  for (const auto& [i, v] : sol) vsol(id_map[i]) = (double) v;
  return getScore(vsol);
}

int Checker::checkLpFeasibility(map<long long, double> &sol){
  VectorXd vsol (n); vsol.fill(0);
  for (const auto& [i, v] : sol) vsol(id_map[i]) = v;
  return checkLpFeasibility(vsol);
}

int Checker::checkIlpFeasibility(map<long long, long long> &sol){
  VectorXd vsol (n); vsol.fill(0);
  for (const auto& [i, v] : sol) vsol(id_map[i]) = (double) v;
  return checkIlpFeasibility(vsol);
}

LsrChecker::~LsrChecker(){
  delete pg;
}

LsrChecker::LsrChecker(const LsrProb &prob, double epsilon): prob(prob), epsilon(epsilon){
  pg = new PgManager();
}

double LsrChecker::getScore(map<long long, double> &sol){
  RMatrixXd out_tuples; vector<long long> out_ids;
  vector<long long> ids; ids.reserve(sol.size());
  for (const auto& [k, v] : sol) ids.push_back(k);
  pg->getSelectedTuples(out_tuples, out_ids, prob.det_sql.table_name, ids, {prob.det_sql.obj_col});
  long long score = 0.0;
  for (int i = 0; i < (int) out_ids.size(); i ++){
    score += out_tuples(0, i)*sol[out_ids[i]];
  }
  if (!prob.det_sql.is_maximize) score = -score;
  return score;
}

double LsrChecker::getScore(map<long long, long long> &sol){
  map<long long, double> dsol;
  for (const auto& [k, v] : sol) dsol[k] = (double) v;
  return getScore(dsol);
}

int LsrChecker::checkLpFeasibility(map<long long, double> &sol){
  RMatrixXd out_tuples; vector<long long> out_ids;
  vector<long long> ids; ids.reserve(sol.size());
  for (const auto& [k, v] : sol){
    if (isLess(v, 0, epsilon)) return LbVariable;
    if (isGreater(v, prob.det_sql.u, epsilon)) return UbVariable;
    ids.push_back(k);
  }
  pg->getSelectedTuples(out_tuples, out_ids, prob.det_sql.table_name, ids, prob.det_sql.att_cols, prob.det_sql.filter_cols, prob.det_sql.filter_intervals);
  int n = (int) out_ids.size();
  // cout << n << " " << ids.size() << endl;
  if (n != (int) ids.size()) return BadFilter;
  int m = (int) prob.det_sql.att_cols.size();
  VectorXd Ax (m); Ax.fill(0);
  double count = 0;
  for (int j = 0; j < m; j ++){
    for (int i = 0; i < n; i ++){
      Ax(j) += out_tuples(j, i)*sol[out_ids[i]];
    }
  }
  for (int i = 0; i < n; i ++){
    count += sol[out_ids[i]];
  }
  for (int i = 0; i < m; i ++){
    if (prob.bl(i) != -DBL_MAX && isLess(Ax(i), prob.bl(i), epsilon)) return LbConstraint;
    if (prob.bu(i) != DBL_MAX && isGreater(Ax(i), prob.bu(i), epsilon)) return UbConstraint;
  }
  if (prob.cl != -DBL_MAX && isLess(count, prob.cl, epsilon)) return LbConstraint;
  if (prob.cu != DBL_MAX && isGreater(count, prob.cu, epsilon)) return UbConstraint;
  return Feasibility;
}

int LsrChecker::checkIlpFeasibility(map<long long, long long> &sol){
  map<long long, double> dsol;
  for (const auto& [k, v] : sol) dsol[k] = (double) v;
  return checkLpFeasibility(dsol);
}