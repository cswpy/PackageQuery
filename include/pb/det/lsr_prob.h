#pragma once

#include "det_bound.h"

#include "pb/det/det_sql.h"
#include "pb/util/udeclare.h"

using namespace pb;

class LsrProb{
public:
  DetSql &det_sql;
  string partition_name;
  VectorXd bl, bu;
  double cl, cu; // Bound on count constraint if exists
  DetBound det_bound;
  int seed;
public:
  ~LsrProb();
  LsrProb(DetSql &det_sql, string partition_name, int seed=-1);
  double generateBounds(double E, double alpha, double hardness);
  void setSeed(int seed);
};