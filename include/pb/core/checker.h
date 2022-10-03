#pragma once

#include "pb/det/det_prob.h"
#include "pb/det/lsr_prob.h"
#include "pb/util/upostgres.h"
#include "pb/util/udeclare.h"

using namespace pb;

static const double kCheckerEpsilon = 1e-6;

class Checker{
private:
  unordered_map<long long, int> id_map;
protected:
  const DetProb &prob;
  int n, m;
  double epsilon;
public:
  Checker(const DetProb &prob, double epsilon=kCheckerEpsilon);
  virtual double getScore(const VectorXd &sol);
  virtual int checkLpFeasibility(const VectorXd &sol);
  virtual int checkIlpFeasibility(const VectorXd &sol);
  double getScore(map<long long, double> &sol);
  double getScore(map<long long, long long> &sol);
  int checkLpFeasibility(map<long long, double> &sol);
  int checkIlpFeasibility(map<long long, long long> &sol);
};

class LsrChecker{
private:
  const LsrProb &prob;
  double epsilon;
  PgManager *pg;
public:
  ~LsrChecker();
  LsrChecker(const LsrProb &prob, double epsilon=kCheckerEpsilon);
  double getScore(map<long long, double> &sol);
  double getScore(map<long long, long long> &sol);
  int checkLpFeasibility(map<long long, double> &sol);
  int checkIlpFeasibility(map<long long, long long> &sol);
};