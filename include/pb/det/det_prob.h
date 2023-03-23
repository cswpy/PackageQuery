#pragma once

#include "det_bound.h"

#include "pb/util/udeclare.h"
#include "pb/det/det_sql.h"

using namespace pb;

// bl <= Ax <= bu
// l <= x <= u
// Maximizing cx
class DetProb{
public:
  DetSql *det_sql;
  RMatrixXd A;
  VectorXd bl, bu, c, l, u;
  vector<long long> ids;
  DetBound det_bound;
  int seed;
public:
  ~DetProb();
  DetProb();
  DetProb(const DetProb& dp);
  void operator=(const DetProb& dp);
  DetProb(int m, int n);
  DetProb(DetSql &det_sql, long long n, int seed=-1);
  void resize(int m, int n);
  void uniformGenerate(int n, int expected_n, double att_var, double outlier_prob, bool restrict_count=true, bool is_positive=false, bool is_translate=false, int seed=-1);
  void normalGenerate(int n, int expected_n, double att_var, double outlier_prob, bool restrict_count=true, int seed=-1);
  double generateBounds(double E, double alpha, double hardness);
  void normalizeObjective();
  void truncate();
  void setSeed(int seed);
  void copyBounds(VectorXd &bl, VectorXd &bu, double cl, double cu);
};