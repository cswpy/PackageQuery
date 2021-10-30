#pragma once

#include <Eigen/Dense>

using namespace Eigen;

class Dual{
public:
  int n, m, iteration_count, mini_iteration_count, status;
  VectorXd x;
  VectorXi bhead;
  double exe_solve, score;
public:
  ~Dual();
  Dual(int core, const MatrixXd& A, const VectorXd& bl, const VectorXd& bu, const VectorXd& c, const VectorXd& l, const VectorXd& u);
};