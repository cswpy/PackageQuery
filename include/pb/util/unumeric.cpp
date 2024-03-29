#include <cstdlib>
#include <boost/math/special_functions/erf.hpp>

#include "unumeric.h"

using boost::math::erf_inv;

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
  return isEqual(frac, 0, eps) || isEqual(frac, 1, eps);
}

long long ceilDiv(long long x, long long q){
  lldiv_t d = lldiv(x, q);
  return d.quot + (d.rem != 0LL);
}

int ceilDiv(int x, int q){
  div_t d = div(x, q);
  return d.quot + (d.rem != 0);
}


double normalQuantile(double u, double v, double p){
  return u + sqrt(2*v) * erf_inv(2*p - 1);
}

double normalCdf(double u, double v, double x){
  return 0.5 + 0.5 * erf((x-u) / (sqrt(2*v)));
}

double pctError(double x, double ground){
  // Assume ground > x
  return (ground - x) / fabs(ground) * 100.0;
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

void MeanVar::add(MeanVar &mv){
  if (attr_count == mv.attr_count){
    if (sample_count){
      long long total_sample_count = sample_count + mv.sample_count;
      VectorXd delta = mv.mean - mean;
      mean = (sample_count * mean + mv.sample_count * mv.mean) / total_sample_count;
      M2 += mv.M2 + delta.cwiseProduct(delta)*sample_count*mv.sample_count/total_sample_count;
      sample_count = total_sample_count;
    } else{
      sample_count = mv.sample_count;
      mean = mv.mean;
      M2 = mv.M2;
    }
  }
}

void MeanVar::reset(){
  mean.fill(0);
  M2.fill(0);
  sample_count = 0;
}

VectorXd MeanVar::getMean(){
  return mean;
}

// Unbiased estimate
VectorXd MeanVar::getVar(){
  if (sample_count <= 1) return M2;
  return M2 / (sample_count - 1);
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

void ScalarMeanVar::add(ScalarMeanVar &smv){
  if (sample_count){
    long long total_sample_count = sample_count + smv.sample_count;
    double delta = smv.mean - mean;
    mean = (sample_count * mean + smv.sample_count * smv.mean) / total_sample_count;
    M2 += smv.M2 + delta*delta*sample_count*smv.sample_count/total_sample_count;
    sample_count = total_sample_count;
  } else{
    // Important to avoid both smv being zero sample_count
    sample_count = smv.sample_count;
    mean = smv.mean;
    M2 = smv.M2;
  }
}

double ScalarMeanVar::getMean(){
  return mean;
}

// Unbiased estimate
double ScalarMeanVar::getVar(){
  if (sample_count <= 1) return 0;
  return M2 / (sample_count - 1);
}

double ScalarMeanVar::getM2(){
  return M2;
}