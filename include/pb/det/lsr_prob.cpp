#include "lsr_prob.h"

#include "pb/util/uconfig.h"
#include "pb/util/upostgres.h"

LsrProb::~LsrProb(){
}

LsrProb::LsrProb(DetSql &det_sql, string partition_name, int seed): det_sql(det_sql), partition_name(partition_name){
  setSeed(seed);
  int m = (int) det_sql.att_cols.size();
  bl.resize(m); bl.fill(-DBL_MAX);
  bu.resize(m); bu.fill(DBL_MAX);
  cl = -DBL_MAX;
  cu = DBL_MAX;
  last_eah = {-1, -1, -1};
  PgManager pg = PgManager();
  // If there is filtering, our means, vars are wrong
  // This leads to infeasible query if the filter ratio is high
  // Hence rebound is needed to correct means and vars for filtering
  Stat* stat = pg.readStats(det_sql.table_name);
  vector<double> means (m);
  vector<double> vars (m);
  for (int i = 0; i < m; i ++){
    int index = stat->getIndex(det_sql.att_cols[i]);
    means[i] = stat->mean(index);
    vars[i] = stat->getVar(index);
  }
  delete stat;
  det_bound = DetBound(det_sql.att_senses, means, vars, LsrProb::seed);
}

double LsrProb::generateBounds(double E, double alpha, double hardness){
  last_eah = {E, alpha, hardness};
  if (det_sql.has_count_constraint){
    cl = E*kLowerCountFactor;
    cu = E*kUpperCountFactor;
  } else {
    cl = -DBL_MAX;
    cu = DBL_MAX;
  }
  bl.fill(-DBL_MAX);
  bu.fill(DBL_MAX);
  return det_bound.sampleHardness(E, alpha, hardness, bl, bu);
}

void LsrProb::setSeed(int seed){
  LsrProb::seed = seed;
  if (seed < 0) {
    random_device rd;
    LsrProb::seed = rd();
  }
  det_bound.setSeed(seed);
}

void LsrProb::rebound(VectorXd means, VectorXd vars){
  det_bound = DetBound(det_sql.att_senses, means, vars, seed);
  auto [E, alpha, hardness] = last_eah;
  if (alpha >= 0) generateBounds(E, alpha, hardness);
}