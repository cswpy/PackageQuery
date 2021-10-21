#pragma once

#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

class Reducer{
public:
  MatrixXd* A;
  VectorXd* b, *c, *u;
  int status;
public:
  ~Reducer();
  Reducer(int core, MatrixXd* AA, VectorXd* bb, VectorXd* cc, VectorXd* uu);
};