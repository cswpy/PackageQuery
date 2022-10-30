#pragma once

#include "pb/util/udeclare.h"
#include "pb/util/unumeric.h"

using namespace pb;

class VG {
public:
  virtual ~VG(){};
  virtual void sample(VectorXd &x, int n, int seed) = 0;
  double sample(int seed=-1){
    VectorXd x; sample(x, 1, seed);
    return x(0);
  }
  double mean(int n, int seed=-1){
    VectorXd x; sample(x, n, seed);
    return x.mean();
  }
  double variance(int n, int seed=-1){
    VectorXd x; sample(x, n, seed);
    ScalarMeanVar smv = ScalarMeanVar();
    for (int i = 0; i < n; i ++) smv.add(x(i));
    return smv.getVar();
  }
};

class LinearVG: public VG {
private:
  vector<VG*> vgs;
  vector<double> ms;
public:
  LinearVG(vector<VG*> vgs, vector<double> ms){
    this->vgs = vgs;
    this->ms = linearCombination(ms);
  }
  void sample(VectorXd &x, int n, int seed){
    x.resize(n);
    vector<VectorXd> multi_samples (vgs.size());
    for (int i = 0; i < (int) vgs.size(); i ++){
      vgs[i]->sample(multi_samples[i], n, seed);
    }
    for (int i = 0; i < n; i ++){
      x(i) = 0;
      for (int j = 0; j < (int) vgs.size(); j ++){
        x(i) += multi_samples[j][i]*ms[j];
      }
    }
  }
};

class MixtureVG: public VG {
private:
  default_random_engine gen;
  vector<VG*> vgs;
  vector<double> ms;
public:
  MixtureVG(vector<VG*> vgs, vector<double> ms){
    this->vgs = vgs;
    this->ms = linearCombination(ms);
  }
  void sample(VectorXd &x, int n, int seed){
    x.resize(n);
    if (seed < 0){
      random_device rd;
      seed = rd();
    }
    gen.seed(seed);
    uniform_real_distribution dist (0.0, 1.0);
    vector<int> counts (vgs.size(), 0);
    for (int i = 0; i < n; i ++){
      double v = dist(gen);
      double s = 0;
      for (int j = 0; j < (int) vgs.size(); j++){
        s += ms[j];
        if (v < s){
          counts[j] ++;
          break;
        }
      }
    }
    int start_index = 0;
    for (int j = 0; j < (int) vgs.size(); j++){
      VectorXd samples;
      vgs[j]->sample(samples, counts[j], seed);
      if (counts[j]) memcpy(&(x(start_index)), &(samples(0)), counts[j]*sizeof(double));
      start_index += counts[j];
    }
    shuffle(x.begin(), x.end(), gen);
  }
};

class EVG: public VG {
protected:
  double vmin, vmax;
public:
  double min(){ return vmin; };
  double max(){ return vmax; };
};

class EEVG: public EVG {
public:
  virtual double quantile(double p) = 0;
};

class Dist: public EEVG{
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