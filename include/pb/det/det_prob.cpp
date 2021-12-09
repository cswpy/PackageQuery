#include "det_prob.h"
#include "pb/util/uconfig.h"
#include "pb/util/udebug.h"

double DetProb::kTranslation = 0;

DetProb::~DetProb(){
}

DetProb::DetProb(){
}

DetProb::DetProb(int m, int n){
  resize(m, n);
}

void DetProb::resize(int m, int n){
  A.resize(m, n);
  bl.resize(m); bl.fill(-DBL_MAX);
  bu.resize(m); bu.fill(DBL_MAX);
  c.resize(n); c.fill(0);
  l.resize(n); l.fill(0);
  u.resize(n); u.fill(1);
}

void DetProb::uniformGenerate(int n, int expected_n, double att_var, double outlier_prob, bool restrict_count, bool is_positive, bool is_translate){
  if (restrict_count) resize(4, n);
  else resize(3, n);
  double left = -sqrt(3*att_var);
  double right = sqrt(3*att_var);
  if (is_positive){
    double translation = 0;
    if (is_translate) translation = kTranslation;
    left = translation * right;
    right = (2 + translation) * right;
  }
  #pragma omp parallel num_threads(kPCore)
  {
    std::random_device rd;
    std::default_random_engine gen (rd());
    std::uniform_real_distribution dist(left, right);
    #pragma omp for
    for (int i = 0; i < n; i ++){
      A(0, i) = dist(gen);
      A(1, i) = dist(gen);
      A(2, i) = dist(gen);
      if (restrict_count) A(3, i) = 1;
      c(i) = dist(gen);
    }
  }
  std::random_device rd;
  std::default_random_engine gen (rd());
  std::normal_distribution dist((left+right)/2*sqrt(expected_n), sqrt(att_var*expected_n));
  //std::uniform_real_distribution dist(left*sqrt(expected_n), right*sqrt(expected_n));
  double tol = sqrt(att_var*expected_n/outlier_prob);
  double v = dist(gen);
  bl(0) = v-tol;
  bu(0) = v+tol;
  bu(1) = dist(gen)+tol;
  bl(2) = dist(gen)-tol;
  if (restrict_count){
    bl(3) = expected_n/2;
    bu(3) = expected_n*3/2;
  }
}

void DetProb::normalGenerate(int n, int expected_n, double att_var, double outlier_prob, bool restrict_count){
  if (restrict_count) resize(4, n);
  else resize(3, n);
  #pragma omp parallel num_threads(kPCore)
  {
    std::random_device rd;
    std::default_random_engine gen (rd());
    std::normal_distribution dist(0.0, sqrt(att_var));
    #pragma omp for
    for (int i = 0; i < n; i ++){
      A(0, i) = dist(gen);
      A(1, i) = dist(gen);
      A(2, i) = dist(gen);
      if (restrict_count) A(3, i) = 1;
      c(i) = dist(gen);
    }
  }
  std::random_device rd;
  std::default_random_engine gen (rd());
  std::normal_distribution dist(0.0, sqrt(expected_n*att_var));
  double tol = sqrt(att_var*expected_n/outlier_prob);
  double v = dist(gen);
  bl(0) = v-tol;
  bu(0) = v+tol;
  bu(1) = dist(gen)+tol;
  bl(2) = dist(gen)-tol;
  if (restrict_count){
    bl(3) = expected_n/2;
    bu(3) = expected_n*3/2;
  }
}