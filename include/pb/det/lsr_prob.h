#pragma once

#include "det_bound.h"

#include "pb/util/udeclare.h"

using namespace pb;

class LsrProb{
public:
  string table_name, partition_name, obj_col;
  bool is_maximize;
  vector<string> cols;
  vector<int> consSense;
  // Assuming upper bounded by u and lower bounded by 0 for all tuples.
  long long u;

  vector<string> filter_cols;
  vector<pair<double, double>> filter_intervals;

  VectorXd bl, bu;
  double cl, cu; // Bound on count constraint if exists
  DetBound detBound;
public:
  ~LsrProb();
  LsrProb();
  LsrProb(string table_name, string partition_name, string obj_col, bool is_maximize, vector<string> cols, vector<int> consSense, long long u=1LL);
  double boundGenerate(double E, double alpha, double hardness);
  void addFilter(string col, double l, double u);
  void setSeed(int seed);
};