#pragma once

#include <queue>
#include <Eigen/Dense>
#include "abstract_walker.h"

using namespace Eigen;
using namespace std;

class PseudoWalker: public AbstractWalker{

public:
  VectorXd p, x;
  vector<double> steps;
  priority_queue<pair<double, int>> pq;
  bool enable_correction;
  int step_count, correction_count;
  double sum_p;

private:
  int executeStep(int i);

public:
  ~PseudoWalker();
  PseudoWalker(VectorXd p, bool enable_correction=true);
  int step();
};
