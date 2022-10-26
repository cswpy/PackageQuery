#include "dist.h"

#include "pb/lib/normal.hpp"
#include "pb/lib/truncated_normal.hpp"
#include "pb/lib/log_normal.hpp"
#include "pb/lib/log_normal_truncated_ab.hpp"

Normal::Normal(double mean, double variance){
  assert(variance > 0);
  u = mean;
  v = variance;
  vmin = -DBL_MAX;
  vmax = DBL_MAX;
  sigma = sqrt(variance);
}

double Normal::pdf(double x){
  return normal_ms_pdf(x, u, sigma);
}

double Normal::cdf(double x){
  return normal_ms_cdf(x, u, sigma);
}

double Normal::quantile(double p){
  if (p <= 0) return vmin;
  if (p >= 1) return vmax;
  return normal_ms_cdf_inv(p, u, sigma);
}

void Normal::sample(VectorXd &x, int n, int seed){
  x.resize(n);
  if (seed < 0){
    random_device rd;
    seed = rd();
  }
  gen.seed(seed);
  normal_distribution<double> distribution(u, sigma);
  for (int i = 0; i < n; i ++){
    x(i) = distribution(gen);
  }
}

double Normal::left_cvar(double p){
  if (p <= 0) return vmin;
  if (p >= 1) return u;
  return truncated_normal_b_mean(u, sigma, quantile(p));
}

////////////////////////////////////////

Uniform::Uniform(double a, double b){
  assert(b > a);
  this->a = a;
  this->b = b;
  u = (a+b)/2;
  v = (b-a)*(b-a)/12.0;
  vmin = a;
  vmax = b;
}

double Uniform::pdf(double x){
  if (x < a || x > b) return 0;
  return 1.0/(b-a);
}

double Uniform::cdf(double x){
  if (x < a) return 0.0;
  if (x > b) return 1.0;
  return (x-a)/(b-a);
}

double Uniform::quantile(double p){
  if (p <= 0) return vmin;
  if (p >= 1) return vmax;
  return a + (b-a)*p;
}

void Uniform::sample(VectorXd &x, int n, int seed){
  x.resize(n);
  if (seed < 0){
    random_device rd;
    seed = rd();
  }
  gen.seed(seed);
  uniform_real_distribution<double> distribution(a, b);
  for (int i = 0; i < n; i ++){
    x(i) = distribution(gen);
  }
}

double Uniform::left_cvar(double p){
  if (p <= 0) return vmin;
  if (p >= 1) return u;
  double tmp = quantile(p);
  return (a + tmp) * (tmp - a) / (2 * p * (b-a));
}

LogNormal::LogNormal(double mu, double sigma){
  assert(sigma > 0);
  this->mu = mu;
  this->sigma = sigma;
  u = log_normal_mean(mu, sigma);
  v = log_normal_variance(mu, sigma);
  vmin = 0;
  vmax = DBL_MAX;
}

double LogNormal::pdf(double x){
  return log_normal_pdf(x, mu, sigma);
}

double LogNormal::cdf(double x){
  return log_normal_cdf(x, mu, sigma);
}

double LogNormal::quantile(double p){
  if (p <= 0) return vmin;
  if (p >= 1) return vmax;
  return log_normal_cdf_inv(p, mu, sigma);
}

void LogNormal::sample(VectorXd &x, int n, int seed){
  x.resize(n);
  if (seed < 0){
    random_device rd;
    seed = rd();
  }
  gen.seed(seed);
  lognormal_distribution<double> distribution(mu, sigma);
  for (int i = 0; i < n; i ++){
    x(i) = distribution(gen);
  }
}

double LogNormal::left_cvar(double p){
  if (p <= 0) return vmin;
  if (p >= 1) return u;
  return log_normal_truncated_ab_mean(mu, sigma, 0, quantile(p));
}