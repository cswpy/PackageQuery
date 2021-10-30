#pragma once

#include <Eigen/Dense>
#include "dual.h"

using namespace Eigen;
using namespace std;

class DualReducer{
public:
  VectorXd best_x;
  int status, layer_count;
  double best_score, exe_solve;
  vector<MatrixXd*> As;
  vector<VectorXd*> bls, bus, cs, ls, us;
  vector<VectorXi*> original_indices;
  vector<Dual*> duals;
public:
  ~DualReducer();
  DualReducer(int core, MatrixXd* AA, VectorXd* bbl, VectorXd* bbu, VectorXd* cc, VectorXd* ll, VectorXd* uu);
};
