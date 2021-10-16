#define FMT_HEADER_ONLY

#include <iostream>
#include <Eigen/Dense>
#include <iostream>
#include <float.h>
#include <bitset>
#include <random>

#include "utility.h"
#include "walker.h"
#include "pseudo_walker.h"
#include "gurobi_lattice_solver.h"
#include "lattice_solver.h"
#include "fmt/core.h"

using namespace std;
using namespace Eigen;

void report(string problem, double time_budget){
  VectorXd sol = readSolution(problem);
  LatticeSolver ls = LatticeSolver(problem);
  fmt::print("Problem:{} n:{} m:{} nnz:{}\n", problem, ls.n, ls.m, ls.lattice_dirs.getSize());
  if (!kIgnoreDegenerate){
    cout << ls.degenerate_count << " degenerate direction out of " << ls.n << endl;
  }
  ls.solve(time_budget);
  fmt::print("Try count:{} Average step count:{:.2Lf} Presolved time:{:.5Lf}ms Relaxed time:{:.5Lf}ms Tableau time:{:.5Lf}ms\n", ls.try_count, ls.avg_step_count, ls.exe_presolve, ls.exe_relaxed, ls.exe_tableau);
  double relaxed_score = ls.getObjValue(ls.uncrushed_r0);
  double obj_score = ls.getObjValue(sol);
  double relax_gap = abs((ls.best_c_score - relaxed_score) / relaxed_score) * 100;
  double obj_gap = abs((ls.best_c_score - obj_score) / obj_score) * 100;
  double lp_gap = abs((obj_score - relaxed_score) / relaxed_score) * 100;
  double l1_dist = (sol-ls.uncrushed_r0).lpNorm<1>();
  bool feas = false;
  if (ls.best_c_score == -DBL_MAX){
    ls.best_c_score = nan("");
    relax_gap = nan("");
    obj_gap = nan("");
  } else feas = ls.checkFeasibility(ls.best_x);
  fmt::print("Best_score:{:.2Lf}\n", ls.best_c_score);
  fmt::print("Relaxed_score:{:.2Lf} gap:{:.2Lf}% L1_dist:{:.2Lf}\n", relaxed_score, relax_gap, l1_dist);
  fmt::print("Obj_score:{:.2Lf} gap:{:.2Lf}% LP_gap:{:.2Lf}%\n", obj_score, obj_gap, lp_gap);
  cout << feas << " " << ls.checkFeasibility(ls.uncrushed_r0) << " " << ls.checkFeasibility(sol) << endl;
}

void report(MatrixXd A, VectorXd b, VectorXd c, double lb_value, double ub_value, double time_budget, bool is_milp){
  LatticeSolver ls = LatticeSolver(A, b, c, lb_value, ub_value);
  fmt::print("n:{} m:{} nnz:{}\n", ls.n, ls.m, ls.lattice_dirs.getSize());
  if (!kIgnoreDegenerate){
    cout << ls.degenerate_count << " degenerate direction out of " << ls.n << endl;
  }
  ls.solve(time_budget);
  fmt::print("Try count:{} Average step count:{:.2Lf} Presolved time:{:.5Lf}ms Relaxed time:{:.5Lf}ms Tableau time:{:.5Lf}ms Init time:{:.5Lf}ms\n", ls.try_count, ls.avg_step_count, ls.exe_presolve, ls.exe_relaxed, ls.exe_tableau, ls.exe_init);
  double relaxed_score = ls.getObjValue(ls.uncrushed_r0);
  double relax_gap = abs((ls.best_c_score - relaxed_score) / relaxed_score) * 100;
  bool feas = false;
  if (ls.best_c_score == -DBL_MAX){
    ls.best_c_score = nan("");
    relax_gap = nan("");
  } else feas = ls.checkFeasibility(ls.best_x);
  fmt::print("Best_score:{:.2Lf}\n", ls.best_c_score);
  fmt::print("Relaxed_score:{:.2Lf} gap:{:.2Lf}%\n", relaxed_score, relax_gap);
  if (is_milp){
    ls.milpSolve();
    double obj_score = ls.milp_score;
    double l1_dist = (ls.x0-ls.uncrushed_r0).lpNorm<1>();
    double obj_gap = abs((ls.best_c_score - obj_score) / obj_score) * 100;
    if (ls.best_c_score == -DBL_MAX) obj_gap = nan("");
    double lp_gap = abs((obj_score - relaxed_score) / relaxed_score) * 100;
    fmt::print("Milp time:{:.5Lf}ms Obj_score:{:.2Lf} gap:{:.2Lf}% LP_gap:{:.2Lf}% L1_dist:{:.2Lf}\n", ls.exe_milp, obj_score, obj_gap, lp_gap, l1_dist);
  }
  cout << feas << " " << ls.checkFeasibility(ls.uncrushed_r0) << endl;
}

void reportGurobi(MatrixXd A, VectorXd b, VectorXd c, double lb_value, double ub_value, double time_budget, bool is_milp){
  GurobiLatticeSolver ls = GurobiLatticeSolver(A, b, c, lb_value, ub_value);
  fmt::print("n:{} m:{} nnz:{}\n", ls.n, ls.m, ls.lattice_dirs.getSize());
  if (!kIgnoreDegenerate){
    cout << ls.degenerate_count << " degenerate direction out of " << ls.n << endl;
  }
  ls.solve(time_budget);
  fmt::print("Try count:{} Average step count:{:.2Lf} Relaxed time:{:.5Lf}ms Tableau time:{:.5Lf}ms Init time:{:.5Lf}ms\n", ls.try_count, ls.avg_step_count, ls.exe_relaxed, ls.exe_tableau, ls.exe_init);
  double relaxed_score = ls.getObjValue(ls.r0);
  double relax_gap = abs((ls.best_c_score - relaxed_score) / relaxed_score) * 100;
  bool feas = false;
  if (ls.best_c_score == -DBL_MAX){
    ls.best_c_score = nan("");
    relax_gap = nan("");
  } else feas = ls.checkFeasibility(ls.best_x);
  fmt::print("Best_score:{:.2Lf}\n", ls.best_c_score);
  fmt::print("Relaxed_score:{:.2Lf} gap:{:.2Lf}%\n", relaxed_score, relax_gap);
  if (is_milp){
    ls.milpSolve();
    double obj_score = ls.milp_score;
    double l1_dist = (ls.x0-ls.r0).lpNorm<1>();
    double obj_gap = abs((ls.best_c_score - obj_score) / obj_score) * 100;
    if (ls.best_c_score == -DBL_MAX) obj_gap = nan("");
    double lp_gap = abs((obj_score - relaxed_score) / relaxed_score) * 100;
    fmt::print("Milp time:{:.5Lf}ms Obj_score:{:.2Lf} gap:{:.2Lf}% LP_gap:{:.2Lf}% L1_dist:{:.2Lf}\n", ls.exe_milp, obj_score, obj_gap, lp_gap, l1_dist);
  }
  cout << feas << " " << ls.checkFeasibility(ls.r0) << endl;
}

void generateProlem(int n, MatrixXd& A, VectorXd& b, VectorXd& c){
  A.resize(6, n); b.resize(6); c.resize(n);
  default_random_engine gen {static_cast<long unsigned int>(time(0))};
  uniform_real_distribution u_dist(0.0, 1.0);
  int expected_numvar = 10;
  normal_distribution n_dist(0.5*expected_numvar, 1.0/12*expected_numvar);
  normal_distribution n_dist_c(0.0, 200.0);
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
  MatrixXd A; VectorXd b, c;
  int n = 1000000;
  generateProlem(n, A, b, c);
  // GurobiLatticeSolver ls = GurobiLatticeSolver(A, b, c, 0, 1);
  // fmt::print("n:{} m:{} nnz:{}\n", ls.n, ls.m, ls.lattice_dirs.getSize());
  // cout << ls.degenerate_count << " degenerate direction out of " << ls.n << endl;
  reportGurobi(A, b, c, 0, 1, 2, true);
  // LatticeSolver ls = LatticeSolver(A, b, c, 0, 1);
  // fmt::print("n:{} m:{} nnz:{}\n", ls.n, ls.m, ls.lattice_dirs.getSize());
  // ls.solve(2);
  // fmt::print("Try count:{} Average step count:{:.2Lf} Presolved time:{:.5Lf}ms Relaxed time:{:.5Lf}ms\n", ls.try_count, ls.avg_step_count, ls.exe_presolve, ls.exe_relaxed);
  // cout << "Sol status:" << solMessage(ls.sol_status) << " Best score:" << ls.best_c_score << endl;
  // cout << ls.best_x.sum() << endl;
  // cout << A << endl;
  // print(b);
  // print(c);
  auto problems = getProblemSizes();
  int cnt = 0;
  for (const auto& pi : problems){
    if (0 <= pi.second && pi.second <= 10000){
      string problem = pi.first;
      cnt ++;
      //GurobiLatticeSolver ls = GurobiLatticeSolver(problem);
      //fmt::print("Problem:{} n:{} m:{} nnz:{}\n", problem, ls.n, ls.m, ls.lattice_dirs.getSize());
      //cout << ls.degenerate_count << " degenerate direction out of " << ls.n << endl;
      //reportGurobi(problem, 10, true);
    }
  }
  // MatrixXd A(5, 3); A << 0, -4, 4, 2, -4, 4, 2, 4, 0, 0, 4, 0, -4, 0, 4;
  // VectorXd b(5); b << 0, 4, 12, 8, 0;
  // VectorXd c(3); c << 0, 0, 1;
  // string problem = "cod105";
  //string problem = "gen-ip054";
  // report(problem, 10, true);
  // LatticeSolver ls = LatticeSolver(A, b, c);
  // ls.solve(2);
  // print(ls.best_x);
  // LatticeSolver ls = LatticeSolver(problem, true);
  // cout << ls.n << " " << ls.m << endl;
  //report(problem, 10, false);
  // int n = 10;
  // vector<double> splits(n-1);
  // default_random_engine gen {static_cast<long unsigned int>(time(0))};
  // uniform_real_distribution u_dist(0.0, 1.0);
  // for (int i = 0; i < n-1; i++) splits[i] = u_dist(gen);
  // VectorXd coefs = bucketSort(splits);
  // print(coefs);
  // for (int i = n-1; i >= 1; i--) coefs(i) -= coefs(i-1);
  // print(coefs);
  // MatrixXd A(5, 3); A << 0, -4, 4, 2, -4, 4, 2, 4, 0, 0, 4, 0, -4, 0, 4;
  // VectorXd b(5); b << 0, 4, 12, 8, 0;
  // VectorXd c(3); c << 1, 0, 0;
  // LatticeSolver ls = LatticeSolver(A, b, c, false);
  // for (int i = 0; i < ls.n; i ++){
  //   cout << ls.degenerate_violation[i].size() << " ";
  // }
  // cout << endl;
  return 0;
}