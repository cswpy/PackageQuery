#pragma once

#include "pb/util/udeclare.h"
#include "pb/util/unumeric.h"

using namespace pb;

const double kRho = 0.5;

class DetBound{
private:
  vector<int> att_senses;
  vector<double> means, vars;
  MatrixXd M;
  int m;
  default_random_engine gen;
public:
  int seed;
  double eps;
private:
  double computeLowerBound(double E, double rho);
  void computeRotationMatrix();
  void assignBound(double rho, const VectorXd& center, VectorXd& bl, VectorXd& bu);
  VectorXd sampleCenter(double E, double rho, double alpha);
public:
  DetBound();
  DetBound(vector<int> att_senses, vector<double> means, vector<double> vars, int seed=-1, double eps=kNumericEps);
  DetBound(vector<int> att_senses, VectorXd means, VectorXd vars, int seed=-1, double eps=kNumericEps);
  double measureHardness(double E, const VectorXd& bl, const VectorXd& bu);
  double minHardness(double& minE, const VectorXd& bl, const VectorXd& bu);
  double sample(double E, double rho, double alpha, VectorXd& bl, VectorXd& bu); // Return hardness
  double sampleHardness(double E, double alpha, double hardness, VectorXd& bl, VectorXd& bu, double rho=kRho); // Return rho
  void setSeed(int seed);
};