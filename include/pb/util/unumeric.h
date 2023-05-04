#pragma once

#include <vector>
#include <utility>
#include <float.h>

#include "Eigen/Dense"

using Eigen::VectorXd;
using std::pair;
using std::vector;

constexpr double kNumericEps = 1e-8;

bool isEqual(double x, double y, double eps=kNumericEps);
bool isLess(double x, double y, double eps=kNumericEps);
bool isGreater(double x, double y, double eps=kNumericEps);
bool isLessEqual(double x, double y, double eps=kNumericEps);
bool isGreaterEqual(double x, double y, double eps=kNumericEps);
bool isInteger(double x, double eps=kNumericEps);

long long ceilDiv(long long x, long long q);
int ceilDiv(int x, int q);

double floorLf(double v, int precision=-1);
double roundLf(double v, int precision=-1);
double ceilLf(double v, int precision=-1);

bool doesIntersect(pair<double, double> a, pair<double, double> b);

vector<double> linearCombination(vector<double> &v);

// Deprecated since normal.hpp is much more optimized
// double normalQuantile(double u, double v, double p);
// double normalCdf(double u, double v, double x);

double pctError(double x, double ground);
double intGap(double x, double ground);

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

// C++ implementation of Welford's running mean-var algorithm
class MeanVar{
private:
  VectorXd mean, M2, max, min;
  int attr_count;
public:
  long long sample_count;
public:
  MeanVar();
  MeanVar(int attr_count);
  void add(const VectorXd &x);
  void add(double *start, int cycle=1);
  void add(MeanVar &mv);
  void reset();
  VectorXd getMean();
  VectorXd getVar();
  VectorXd getM2();
  VectorXd getRange();
  VectorXd getMin();
  VectorXd getMax();
};

class ScalarMeanVar{
private:
  double mean, M2, min, max;
public:
  long long sample_count;
public:
  ScalarMeanVar();
  ScalarMeanVar(const ScalarMeanVar &smv);
  void add(double x);
  void add(ScalarMeanVar &smv);
  void reset();
  double getMean();
  double getVar();
  double getM2();
  double getRange();
  double getMin();
  double getMax();
  double getBiasedVar();
};