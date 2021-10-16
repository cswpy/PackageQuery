#define FMT_HEADER_ONLY

#include <iostream>
#include <Eigen/Dense>
#include <iostream>
#include <float.h>
#include <random>

#include "utility.h"
#include "lattice_solver.h"
#include "gurobi_lattice_solver.h"
#include "fmt/core.h"

using namespace std;
using namespace Eigen;

#define ta(row, col) (ls.simplex->tableau[ls.simplex->numcols*(row)+(col)])

void generateProlem(int n, MatrixXd& A, VectorXd& b, VectorXd& c){
  A.resize(6, n); b.resize(6); c.resize(n);
  //default_random_engine gen;
  default_random_engine gen {static_cast<long unsigned int>(time(0))};
  uniform_real_distribution u_dist(0.0, 1.0);
  int expected_numvar = 10;
  normal_distribution n_dist(0.5*expected_numvar, 1.0/12*expected_numvar);
  normal_distribution n_dist_c(0.0, 200.0);
  #pragma omp parallel for num_threads(CORE_COUNT)
  for (int i = 0; i < n; i ++){
    A(0, i) = u_dist(gen);
    A(1, i) = -A(0, i);
    A(2, i) = u_dist(gen);
    A(3, i) = -A(2, i);
    A(4, i) = 1;
    A(5, i) = -1;
    c(i) = n_dist_c(gen);
  }
  double tol_scale = 5;
  double tol1 = u_dist(gen) * tol_scale / 2;
  double tol2 = u_dist(gen) * tol_scale / 2;
  b(0) = n_dist(gen) + tol1;
  b(1) = tol1 - b(0);
  b(2) = n_dist(gen) + tol2;
  b(3) = tol2 - b(2);
  b(4) = 20;
  b(5) = -5;
}

int main(){
  int n = 500000;
  int m = 6;
  MatrixXd A (m, n);
  VectorXd b (m);
  VectorXd c (n);
  generateProlem(n, A, b, c);
  VectorXd u (n); u.fill(1);
  LatticeSolver ls = LatticeSolver(80, A, b, c, u);
  // ls.lattice_dirs.print();
  // print(ls.fcv->sample());
  // ls.solve(2);
  ////////////////////
  // MatrixXd A(5, 3); A << 0, -4, 4, 2, -4, 4, 2, 4, 0, 0, 4, 0, -4, 0, 4;
  // VectorXd b(5); b << 0, 4, 12, 8, 0;
  // VectorXd c(3); c << 1, 0, 0;
  // VectorXd u (3); u.fill(10);
  // LatticeSolver ls = LatticeSolver(80, A, b, c, u);
  // ls.lattice_dirs.print();
  GurobiLatticeSolver gls = GurobiLatticeSolver(A, b, c, 0, 1);
  double error = 0;
  for (int i = 0; i < m; i ++){
    int basic_col = gls.bhead[i];
    int my_row = ls.inv_bhead[basic_col];
    for (int j = 0; j < n+m; j ++){
      error += abs(ta(my_row+2, basic_col)-gls.tableau(i, basic_col));
    }
  }
  error /= (m*(m+n));
  cout << error << endl;
  // gls.solve(2);
  // if (gls.status == LS_FOUND) ls.compareReport(gls.best_x);
  // else ls.report();
  // cout << gls.best_c_score << " " << solMessage(gls.status) << endl;
}