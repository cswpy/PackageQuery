#pragma once

#include "pb/util/udeclare.h"

using namespace pb;

// bl <= Ax <= bu
// l <= x <= u
// Maximizing cx
class DetProb{
public:
  static double kTranslation;
  RMatrixXd A;
  VectorXd bl, bu, c, l, u;
public:
  ~DetProb();
  DetProb();
  DetProb(int m, int n);
  void resize(int m, int n);
  void uniformGenerate(int n, int expected_n, double att_var, double outlier_prob, bool restrict_count=true, bool is_positive=false, bool is_translate=false);
  void normalGenerate(int n, int expected_n, double att_var, double outlier_prob, bool restrict_count=true);
};