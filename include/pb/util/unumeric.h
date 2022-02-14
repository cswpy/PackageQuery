#pragma once

constexpr double kNumericEps = 1e-12;

bool isEqual(double x, double y, double eps=kNumericEps);
bool isLess(double x, double y, double eps=kNumericEps);
bool isGreater(double x, double y, double eps=kNumericEps);
bool isLessEqual(double x, double y, double eps=kNumericEps);
bool isGreaterEqual(double x, double y, double eps=kNumericEps);
bool isInteger(double x, double eps=kNumericEps);

long long ceilDiv(long long x, long long q);

template <typename T>
double sign(T value) {
  if (value < 0) return -1.0;
  if (value > 0) return 1.0;
  return 0.0;
}

template <typename T>
double nonNegativeSign(T value) {
  if (value < 0) return -1.0;
  return 1.0;
}

#include "Eigen/Dense"

using Eigen::VectorXd;

// C++ implementation of Welford's running mean-var algorithm
class MeanVar{
private:
  VectorXd mean, M2;
  int attr_count;
public:
  int sample_count;
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
public:
  int sample_count;
public:
  ScalarMeanVar();
  void add(double x);
  void reset();
  double getMean();
  double getVar();
  double getM2();
};