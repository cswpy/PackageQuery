#pragma once

#include "pb/det/det_prob.h"
#include "pb/util/udeclare.h"
#include "boost/dynamic_bitset.hpp"

using namespace pb;

class Dual{
public:
  int n, m, iteration_count, mini_iteration_count, status;
  VectorXd sol;
  VectorXi bhead;
  MatrixXd Binv;
  boost::dynamic_bitset<> inv_bhead;
  double exe_solve, score;
private:
  // Internal members
  VectorXd bu, bl, beta, d;
  VectorXd rho_r, alpha_r, alpha_q, tau;
public:
  ~Dual();
  Dual(int core, const DetProb &prob);
};