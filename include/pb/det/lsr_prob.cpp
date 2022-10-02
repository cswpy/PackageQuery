#include "lsr_prob.h"

#include "pb/util/uconfig.h"
#include "pb/util/upostgres.h"

LsrProb::~LsrProb(){
}

LsrProb::LsrProb(){
}

LsrProb::LsrProb(string table_name, string partition_name, string obj_col, bool is_maximize, vector<string> cols, vector<int> consSense, long long u)
  :table_name(table_name), partition_name(partition_name), obj_col(obj_col), is_maximize(is_maximize), cols(cols), consSense(consSense), u(u){
  int m = (int) cols.size();
  bl.resize(m); bl.fill(-DBL_MAX);
  bu.resize(m); bu.fill(DBL_MAX);
  cl = -DBL_MAX;
  cu = DBL_MAX;
  PgManager pg = PgManager();
  Stat* stat = pg.readStats(table_name);
  vector<double> means (m);
  vector<double> vars (m);
  for (int i = 0; i < m; i ++){
    int index = stat->getIndex(cols[i]);
    means[i] = stat->mean(index);
    vars[i] = stat->getVar(index);
  }
  delete stat;
  detBound = DetBound(consSense, means, vars, kGlobalSeed);
}

double LsrProb::boundGenerate(double E, double alpha, double hardness){
  cl = E/2.0;
  cu = E*3/2.0;
  return detBound.sampleHardness(E, alpha, hardness, bl, bu);
}

void LsrProb::addFilter(string col, double l, double u){
  if (l == -DBL_MAX && u == DBL_MAX) return;
  filter_cols.push_back(col);
  filter_intervals.emplace_back(l, u);
}

void LsrProb::setSeed(int seed){
  detBound.setSeed(seed);
}