#pragma once

#include "pb/det/det_prob.h"
#include "Eigen/Dense"

using Eigen::VectorXd;

class Checker{
protected:
  const DetProb &prob;
  int n, m;
  double epsilon;
public:
  Checker(const DetProb &prob, double epsilon=1e-10);
  virtual double getScore(const VectorXd &sol);
  virtual int checkLpFeasibility(const VectorXd &sol);
  virtual int checkIlpFeasibility(const VectorXd &sol);
};