#define FMT_HEADER_ONLY

#include <iostream>
#include <Eigen/Dense>
#include <iostream>
#include <float.h>
#include <queue>
#include <random>
#include <chrono>
#include <numeric>

#include "map_sort.h"
#include "utility.h"
#include "parallel_pq.h"
#include "fmt/core.h"
#include "fitsio.h"
#include "omp.h"
#include "lattice_solver.h"
#include "pseudo_walker.h"
#include "simplex.h"
#include "gurobi_lattice_solver.h"
#include "gurobi_solver.h"
#include "dual.h"
#include "dual_reducer.h"
#include "reducer.h"
#include "fmt/core.h"

using namespace std;
using namespace Eigen;

void generateQuery(int qi, int n, RMatrixXd& A, VectorXd& bl, VectorXd& bu, VectorXd& c){
  double outlier_prob = 0.6;
  vector<double> lbs = {-1, -1, -1, -1, 0, 0, 0, 0};
  vector<int> expected_ns = {20, 20, 500, 500, 20, 20, 500, 500};
  vector<double> att_vars = {1, 100, 1, 100, 1, 100, 1, 100};
  int expected_n = expected_ns[qi];
  double att_var = att_vars[qi];
  A.resize(4, n); bl.resize(4); bu.resize(4); c.resize(n); 
  default_random_engine gen {static_cast<long unsigned int>(time(0))};
  double left = lbs[qi];
  double right = 1.0;
  uniform_real_distribution u_dist(left, right);
  double multiplicity = att_var / (right-left) * 12;
  //normal_distribution att_dist(10.0, att_var);
  int expected_numvar = expected_n;
  double mean = 0.5*multiplicity*expected_numvar;
  double var = att_var*expected_numvar;
  normal_distribution n_dist(0.0, var);
  //normal_distribution n_dist_c(0.0, 1.0/12);
  uniform_real_distribution n_dist_c(left, right);
  #pragma omp parallel for num_threads(CORE_COUNT)
  for (int i = 0; i < n; i ++){
    A(0, i) = u_dist(gen)*multiplicity;
    A(1, i) = u_dist(gen)*multiplicity;
    A(2, i) = u_dist(gen)*multiplicity;
    A(3, i) = 1;
    c(i) = n_dist_c(gen)*multiplicity;
    //c(i) = A(0, i) + A(1, i) + A(2, i);
  }
  double tol = var*sqrt(1/outlier_prob);
  bl(0) = mean - tol;
  bu(0) = mean + tol;
  bl(1) = -DBL_MAX;
  bu(1) = mean + tol;
  bl(2) = mean- tol;
  bu(2) = DBL_MAX;
  bl(3) = expected_n*1.0/2;
  bu(3) = expected_n*3.0/2;
}

void testAccuracy(){
  int T = 3;
  vector<int> szs = {1000, 10000, 100000};
  vector<int> repeats = {100, 50, 25};
  RMatrixXd A;
  VectorXd bl, bu, c;
  VectorXd u;
  VectorXd l;
  for (int qi = 0; qi < 8; qi ++){
    for (int i = 0; i < T; i ++){
      if (i != 0) continue;
      int fail0 = 0;
      int fail1 = 0;
      double obj_gap0 = 0;
      double obj_gap1 = 0;
      double ap0 = 0;
      double ap1 = 0;
      int succ = 0;
      int n = szs[i];
      for (int r = 0; r < repeats[i]; r++){
        generateQuery(qi, n, A, bl, bu, c);
        u.resize(n); u.fill(1);
        l.resize(n); l.fill(0);
        DualReducer dr0 = DualReducer(16, &A, &bl, &bu, &c, &l, &u, 0);
        DualReducer dr1 = DualReducer(16, &A, &bl, &bu, &c, &l, &u, 1);
        double base_score = dr0.duals[0]->score;
        if (dr0.status != LS_FOUND || dr1.status != LS_FOUND){
          if (dr0.status != LS_FOUND) fail0 ++;
          if (dr1.status != LS_FOUND) fail1 ++;
        } else{
          succ ++;
          double gap0 = (base_score - dr0.best_score) / base_score;
          double gap1 = (base_score - dr1.best_score) / base_score;
          obj_gap0 += gap0;
          obj_gap1 += gap1;
          ap0 += base_score / dr0.best_score;
          ap1 += base_score / dr1.best_score;
        }
      }
      fmt::print("{},{},{:.2Lf},{:.2Lf},{:.10Lf},{:.10Lf},{:.10Lf},{:.10Lf}\n", qi, n, fail0/(double)repeats[i], fail1/(double)repeats[i], obj_gap0/succ, obj_gap1/succ, ap0/succ, ap1/succ);
      //fmt::print("At ({},{},{}) fail_prob({:.2Lf},{:.2Lf}) avg_obj_gap({:.10Lf},{:.10Lf}) avg_ap({:.10Lf},{:.10Lf})\n", qi, i, r, fail0/(double)repeats[i], fail1/(double)repeats[i], obj_gap0/succ, obj_gap1/succ, ap0/succ, ap1/succ);
    }
  }
}

void testDual(){
  int T = 5;
  vector<int> szs = {1000, 10000, 100000, 1000000, 10000000};
  RMatrixXd A;
  VectorXd bl, bu, c;
  VectorXd u;
  VectorXd l;
  for (int qi = 0; qi < 8; qi ++){
    if (qi != 4) continue;
    for (int i = 0; i < T; i ++){
      int n = szs[i];
      generateQuery(qi, n, A, bl, bu, c);
      u.resize(n); u.fill(1);
      l.resize(n); l.fill(0);
      Dual dual = Dual(8, A, bl, bu, c, l, u);
      GurobiSolver gs = GurobiSolver(A, bl, bu, c, l, u);
      gs.solveRelaxed();
      fmt::print("{},{},{:.2Lf},{:.2Lf}\n", qi, n, dual.exe_solve, gs.exe_relaxed+gs.exe_init);
    }
  }
}

void testReducer(){
  int T = 4;
  vector<int> szs = {1000, 10000, 100000, 1000000};
  RMatrixXd A;
  VectorXd bl, bu, c;
  VectorXd u;
  VectorXd l;
  for (int qi = 0; qi < 8; qi ++){
    for (int i = 0; i < T; i ++){
      if (i != 1) continue;
      int n = szs[i];
      generateQuery(qi, n, A, bl, bu, c);
      u.resize(n); u.fill(1);
      l.resize(n); l.fill(0);
      DualReducer dr = DualReducer(8, &A, &bl, &bu, &c, &l, &u, 0);
      double r0_obj = dr.duals[0]->score;
      double gap = (r0_obj - dr.best_score) / r0_obj;
      fmt::print("{},{},{:.2Lf},{:.2Lf},{:.5Lf}\n", qi, n, dr.exe_solve, dr.duals[0]->exe_solve, gap);
    }
  }
}

void testMajorMatrix(){
  vector<string> names = {"0", "1", "2", "3", "4", "5"};
  Profiler pro = Profiler(names);
  int m = 25;
  int N = 10000000;
  Matrix<double, Dynamic, Dynamic, ColMajor> mat1 (m, N); mat1.fill(2);
  RMatrixXd mat2 (m, N); mat2.fill(2);
  RMatrixXd mat3 (m, N); mat3.fill(2);
  int core = 8;
  #pragma omp parallel num_threads(core)
  {
    pro.clock(0);
    for (int i = 0; i < m; i ++){
      #pragma omp for nowait
      for (int j = 0; j < N; j ++){
        mat1(i, j) *= 2;
      }
    }
    #pragma omp barrier
    pro.stop(0);
    pro.clock(1);
    for (int i = 0; i < m; i ++){
      #pragma omp for nowait
      for (int j = 0; j < N; j ++){
        mat2(i, j) *= 2;
      }
    }
    #pragma omp barrier
    pro.stop(1);
    pro.clock(2);
    #pragma omp for nowait
    for (int j = 0; j < N; j ++){
      for (int i = 0; i < m; i ++){
        mat1(i, j) *= 2;
      }
    }
    #pragma omp barrier
    pro.stop(2);
    pro.clock(3);
    #pragma omp for nowait
    for (int j = 0; j < N; j ++){
      for (int i = 0; i < m; i ++){
        mat2(i, j) *= 2;
      }
    }
    #pragma omp barrier
    pro.stop(3);
  }
  pro.clock(4);
  mat3 *= 2;
  pro.stop(4);
  // pro.clock(5);
  // for (int i = 0; i < m; i ++){
  //   for (int j = 0; j < N; j ++){
  //     mat3(i, j) *= 2;
  //   }
  // }
  // pro.stop(5);
  pro.print();
}

int main(){
  testAccuracy();
  //testDual();
  //testReducer();
  //testMajorMatrix();
}