#pragma once

#include "pb/core/vg.h"
#include "pb/util/udeclare.h"

using namespace pb;

class Normal : public Dist{
private:
  double sigma;
public:
  Normal(double mean, double variance);
  double pdf(double x);
  double cdf(double x);
  double quantile(double p);
  void sample(VectorXd &x, int n, int seed=-1);
  double left_cvar(double p);
};

class Uniform : public Dist{
private:
  double a, b;
public:
  Uniform(double a, double b);
  double pdf(double x);
  double cdf(double x);
  double quantile(double p);
  void sample(VectorXd &x, int n, int seed=-1);
  double left_cvar(double p);
};

class LogNormal : public Dist{
private:
  double mu, sigma;
public:
  LogNormal(double mu, double sigma);
  double pdf(double x);
  double cdf(double x);
  double quantile(double p);
  void sample(VectorXd &x, int n, int seed=-1);
  double left_cvar(double p);
};