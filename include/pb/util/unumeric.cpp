#include <cstdlib>

#include "unumeric.h"

bool isEqual(double x, double y, double eps){
  return fabs(x-y) < eps;
}

bool isLessEqual(double x, double y, double eps){
  return x-y < eps;
}

bool isGreaterEqual(double x, double y, double eps){
  return x-y > -eps;
}

bool isLess(double x, double y, double eps){
  return x-y < -eps;
}

bool isGreater(double x, double y, double eps){
  return x-y > eps;
}

bool isInteger(double x, double eps){
  double whole, frac;
  frac = modf(x, &whole);
  return isEqual(frac, 0, eps);
}

long long ceilDiv(long long x, long long q){
  lldiv_t d = lldiv(x, q);
  return d.quot + (d.rem != 0LL);
}

MeanVar::MeanVar(){
}


MeanVar::MeanVar(int attr_count){
  mean.resize(attr_count); mean.fill(0);
  M2.resize(attr_count); M2.fill(0);
  this->attr_count = attr_count;
  sample_count = 0;
}

void MeanVar::add(const VectorXd& x){
  sample_count ++;
  VectorXd delta = x - mean;
  mean += delta / sample_count;
  M2 += delta.cwiseProduct(x - mean);
}

void MeanVar::add(double *start, int cycle){
  sample_count ++;
  for (int i = 0; i < attr_count; i ++){
    double x = start[i*cycle];
    double delta = x - mean(i);
    mean(i) += delta / sample_count;
    M2(i) += delta * (x - mean(i));
  }
}

VectorXd MeanVar::getMean(){
  return mean;
}

VectorXd MeanVar::getVar(){
  if (sample_count == 0) return M2;
  return M2 / sample_count;
}

VectorXd MeanVar::getM2(){
  return M2;
}

ScalarMeanVar::ScalarMeanVar(){
  mean = 0;
  M2 = 0;
  sample_count = 0;
}

void ScalarMeanVar::reset(){
  mean = 0;
  M2 = 0;
  sample_count = 0;
}

void ScalarMeanVar::add(double x){
  sample_count ++;
  double delta = x - mean;
  mean += delta / sample_count;
  M2 += delta * (x - mean);
}

double ScalarMeanVar::getMean(){
  return mean;
}

double ScalarMeanVar::getVar(){
  if (sample_count == 0) return 0;
  return M2 / sample_count;
}

double ScalarMeanVar::getM2(){
  return M2;
}