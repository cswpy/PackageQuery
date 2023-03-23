#include <cstdlib>
#include <boost/math/special_functions/erf.hpp>

#include "unumeric.h"
#include "fmt/core.h"

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

double floorLf(double v, int precision){
  if (precision < 0) return v;
  double scale = pow(10.0, (double) precision);
  return floor(v*scale) / scale;
}

double roundLf(double v, int precision){
  if (precision < 0) return v;
  double scale = pow(10.0, (double) precision);
  return round(v*scale) / scale;
}

double ceilLf(double v, int precision){
  if (precision < 0) return v;
  double scale = pow(10.0, (double) precision);
  return ceil(v*scale) / scale;
}

bool doesIntersect(pair<double, double> a, pair<double, double> b){
  return !(a.first >= b.second || b.first >= a.second);
}

vector<double> linearCombination(vector<double> &v){
  double sum = 0;
  vector<double> res (v.size());
  for (auto x : v){
    assert(x >= 0);
    sum += x;
  }
  assert(sum > 0);
  for (int i = 0; i < (int) v.size(); i ++){
    res[i] = v[i] / sum;
  }
  return res;
}


// double normalQuantile(double u, double v, double p){
//   return u + sqrt(2*v) * erf_inv(2*p - 1);
// }

// double normalCdf(double u, double v, double x){
//   return 0.5 + 0.5 * erf((x-u) / (sqrt(2*v)));
// }

double pctError(double x, double ground){
  // Assume ground > x
  return (ground - x) / fabs(ground) * 100.0;
}

double intGap(double x, double ground){
  // Assume x, ground have same sign
  double fx = fabs(x)+1e-1;
  double fg = fabs(ground)+1e-1;
  return std::max(fx/fg, fg/fx);
}

MeanVar::MeanVar(){
}


MeanVar::MeanVar(int attr_count){
  mean.resize(attr_count); mean.fill(0);
  M2.resize(attr_count); M2.fill(0);
  max.resize(attr_count); max.fill(-DBL_MAX);
  min.resize(attr_count); min.fill(DBL_MAX);
  this->attr_count = attr_count;
  sample_count = 0;
}

void MeanVar::add(const VectorXd& x){
  sample_count ++;
  VectorXd delta = x - mean;
  mean += delta / sample_count;
  M2 += delta.cwiseProduct(x - mean);
  max = max.cwiseMax(x);
  min = min.cwiseMin(x);
}
// max min are not implemented on this method
void MeanVar::add(double *start, int cycle){
  sample_count ++;
  for (int i = 0; i < attr_count; i++){
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
      max = max.cwiseMax(mv.max);
      min = min.cwiseMin(mv.min);
    } else{
      sample_count = mv.sample_count;
      mean = mv.mean;
      M2 = mv.M2;
      max = mv.max;
      min = mv.min;
    }
  }
}

void MeanVar::reset(){
  mean.fill(0);
  M2.fill(0);
  sample_count = 0;
  max.fill(0);
  min.fill(0);
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

VectorXd MeanVar::getRange(){
  return max-min;
}

VectorXd MeanVar::getMin(){
  return min;
}

VectorXd MeanVar::getMax(){
  return max;
}

ScalarMeanVar::ScalarMeanVar(){
  mean = 0;
  M2 = 0;
  sample_count = 0;
  min = DBL_MAX;
  max = -DBL_MAX;
}

ScalarMeanVar::ScalarMeanVar(const ScalarMeanVar &smv){
  mean = smv.mean;
  M2 = smv.M2;
  sample_count = smv.sample_count;
  min = smv.min;
  max = smv.max;
}

void ScalarMeanVar::reset(){
  mean = 0;
  M2 = 0;
  sample_count = 0;
  max = -DBL_MAX;
  min = DBL_MAX;
}

void ScalarMeanVar::add(double x){
  sample_count ++;
  double delta = x - mean;
  mean += delta / sample_count;
  M2 += delta * (x - mean);
  if (x<min) min = x;
  if (x>max) max = x;
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

double ScalarMeanVar::getRange(){
  return max-min;
}

double ScalarMeanVar::getBiasedVar(){
  if (sample_count <= 1) return 0;
  return M2 / sample_count;
}

double ScalarMeanVar::getMin(){
  return min;
}

double ScalarMeanVar::getMax(){
  return max;
}