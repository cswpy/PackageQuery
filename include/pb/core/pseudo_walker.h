#pragma once

#include <queue>
#include "pb/util/udeclare.h"
#include "parallel_pq.h"

using namespace pb;
using std::priority_queue;

class PseudoWalker{

private:
  VectorXd p, x, steps;
  priority_queue<pair<double, int>> pq;
  ParallelPQ* ppq;
  bool enable_correction;
  int n, step_count, correction_count, core;
  double sum_p;

private:
  int executeStep(int i);

public:
  ~PseudoWalker();
  PseudoWalker(VectorXd& p, bool enable_correction=true, int core=1);
  int step();
};
