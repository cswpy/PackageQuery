#pragma once

#include "Eigen/Dense"

using Eigen::VectorXd;

// C++ implementation of Welford's running mean-var algorithm
class MeanVar{
private:
  VectorXd mean, M2;
  int sample_count, attr_count;
public:
  MeanVar();
  MeanVar(int attr_count);
  void add(const VectorXd& x);
  void add(double *start, int cycle=1);
  VectorXd getMean();
  VectorXd getVar();
  VectorXd getM2();
};

class ScalarMeanVar{
private:
  double mean, M2;
  int sample_count;
public:
  ScalarMeanVar();
  void add(double x);
  void reset();
  double getMean();
  double getVar();
  double getM2();
};