#include <iostream>


#include "pseudo_walker.h"
#include "utility.h"

using namespace std;

PseudoWalker::~PseudoWalker(){
}

PseudoWalker::PseudoWalker(VectorXd p, bool enable_correction){
  int dim = (int) p.size();
  assert(dim > 1);
  this->p = p / p.norm();
  double S = this->p.cwiseAbs().sum();
  vector<pair<double, int>> init(dim);
  steps.resize(dim);
  for (int i = 0; i < dim; i++) {
    double abs_pi = abs(this->p(i));
    if (abs_pi > kFloatEps) {
      steps[i] = 2 - 2 * abs_pi * (abs_pi - (S - abs_pi) / (dim - 1.0));
      init[i] = { this->p(i) * this->p(i), i };
    }
  }
  pq = priority_queue(init.begin(), init.end());
  this->enable_correction = enable_correction;
  step_count = 0;
  correction_count = 0;
  sum_p = 0;
  x.resize(dim); x.fill(0);
}

int PseudoWalker::executeStep(int i){
  step_count ++;
  double s = sign(p(i));
  x(i) += s;
  sum_p += abs(p(i));
  if (enable_correction && step_count >= p.size() / log2(p.size()) * (correction_count+1)){
    int dim = p.size();
    vector<pair<double, int>> init(dim);
    for (int i = 0; i < dim; i++) {
      double abs_pi = abs(this->p(i));
      if (abs_pi > kFloatEps) {
        init[i] = { abs_pi * (abs_pi + 2*sum_p) - 2*abs(x(i)), i };
      }
    }
    pq = priority_queue(init.begin(), init.end());
    correction_count ++;
  }
  return s*(i+1);
}

int PseudoWalker::step() {
  const auto& pi = pq.top();
  int i = pi.second;
  double new_pv = pi.first - steps[i];
  pq.pop();
  pq.emplace(new_pv, i);
  return executeStep(i);
}