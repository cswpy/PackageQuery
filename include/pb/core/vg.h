#pragma once

#include "pb/util/udeclare.h"

using namespace pb;

class VG {
public:
  virtual void sample(VectorXd &x, int n, int seed) = 0;
};

class EVG: public VG {
protected:
  double vmin, vmax;
public:
  double min(){ return vmin; };
  double max(){ return vmax; };
  virtual double quantile(double p) = 0;
};

class Dist: public EVG{
protected:
  default_random_engine gen;
  double u, v;
public:
  double mean(){ return u;};
  double variance(){ return v;};
  virtual double pdf(double x) = 0;
  virtual double cdf(double x) = 0;
  virtual double left_cvar(double p) = 0;
  double right_cvar(double p){
    if (p <= 0) return vmax;
    if (p >= 1) return u;
    return (u - left_cvar(1-p)*(1-p)) / p; 
  }
};