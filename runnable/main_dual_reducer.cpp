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

void generateBoundedProlem(int expected_n, double outlier_prob, double att_var, int n, MatrixXd& A, VectorXd& bl, VectorXd& bu, VectorXd& c){
  A.resize(4, n); bl.resize(4); bu.resize(4); c.resize(n); 
  default_random_engine gen {static_cast<long unsigned int>(time(0))};
  // uniform_real_distribution u_dist(0.0, 1.0);
  normal_distribution att_dist(0.0, att_var);
  int expected_numvar = expected_n;
  double mean = 0;
  double var = att_var*expected_numvar;
  normal_distribution n_dist(0.0, var);
  normal_distribution n_dist_c(0.0, 200.0);
  #pragma omp parallel for num_threads(CORE_COUNT)
  for (int i = 0; i < n; i ++){
    A(0, i) = att_dist(gen);
    A(1, i) = att_dist(gen);
    A(2, i) = att_dist(gen);
    A(3, i) = 1;
    c(i) = n_dist_c(gen);
  }
  double tol = var*sqrt(1/outlier_prob);
  bl(0) = mean - tol;
  bu(0) = mean + tol;
  bl(1) = -DBL_MAX;
  bu(1) = mean + tol;
  bl(2) = mean- tol;
  bu(2) = DBL_MAX;
  bl(3) = expected_n/2;
  bu(3) = expected_n*2;
}

void centroid_test(){
  int core = 16;
  int n = 2;
  int m = 4;
  MatrixXd A (m, n); A << 1.0/3, 1, 4.0/3, -1, -1.0/4, 1, 1, -1;
  VectorXd bl (m); bl << 8.0/3, -DBL_MAX, 1.0/2, 0;
  VectorXd bu (m); bu << DBL_MAX, 17.0/3, 3, DBL_MAX;
  VectorXd c (n); c << -1, 2;
  VectorXd l (n); l << 0, 0;
  VectorXd u (n); u << 10, 10;

  Dual* dual = new Dual(core, A, bl, bu, c, l, u);
  print(dual->x);

  if (dual->status == LS_FOUND){
    VectorXd centroid_dir (n); centroid_dir.fill(0);
    VectorXd r0 (n);
    VectorXd stay_r0 (n);
    VectorXd sum_row(m); sum_row.fill(0);
    int stay_count = 0;
    #pragma omp parallel num_threads(core)
    {
      double whole, frac;
      int local_stay_count = 0;
      #pragma omp for nowait
      for (int i = 0; i < n; i ++){
        r0(i) = dual->x(i);
        frac = modf(dual->x(i), &whole);
        if (isEqual(frac, 0)) stay_r0(i) = whole;
        else{
          stay_r0(i) = -1;
          local_stay_count ++;
        }
      }
      #pragma omp atomic
      stay_count += local_stay_count;
      #pragma omp barrier
      for (int j = 0; j < m; j ++){
        double local_sum_row = 0;
        #pragma omp for nowait
        for (int i = 0; i < m+n; i ++){
          if (!dual->inv_bhead[i]){
            // Non-basic
            if (i < n){
              if (isEqual(dual->x(i), dual->l(i))) local_sum_row -= dual->A(j, i);
              else local_sum_row += dual->A(j, i);
            } else{
              int index = i - n;
              if (index == j){
                if (isEqual(dual->x(i), dual->bl(index))) local_sum_row --;
                else local_sum_row ++;
              }
            }
          }
        }
        #pragma omp atomic
        sum_row(j) += local_sum_row;
      }
      #pragma omp for nowait
      for (int i = 0; i < n; i ++){
        if (!dual->inv_bhead[i]){
          // Non-basic
          if (isEqual(dual->x(i), dual->l(i))) centroid_dir(i) = 1.0/n;
          else centroid_dir(i) = -1.0/n;
        }
      }
      #pragma omp barrier
      #pragma omp master
      {
        print(sum_row);
        for (int i = 0; i < m; i ++){
          if (dual->bhead(i) < n){
            double val = 0;
            for (int j = 0; j < m; j ++){
              val += dual->Binv(i, j) * sum_row(j);
            }
            centroid_dir(dual->bhead(i)) = val / n;
          }
        }
      }
      #pragma omp barrier
      #pragma omp single
      {
        print(centroid_dir);
      }
    }
  }
}

void quickRun(){
  int n = 100;
  MatrixXd A;
  VectorXd bl, bu, c;
  generateBoundedProlem(10, 0.8, 3.0, n, A, bl, bu, c);
  VectorXd u (n); u.fill(1);
  VectorXd l (n); l.fill(0);

  // VectorXd* a = new VectorXd(10);
  // a->fill(10);
  // cout << (*a)(1) << endl;
  // VectorXd* b = new VectorXd(10);
  // memcpy(&(*b)(0), &(*a)(0), 10*8);
  // (*a)(1) = 12;
  // print(*a);
  // print(*b);

  GurobiSolver gs = GurobiSolver(A, bl, bu, c, l, u);
  gs.solveRelaxed();
  // print(gs.r0);
  cout << gs.exe_relaxed << " " << gs.iteration_count << endl;
  cout << gs.relaxed_cscore << endl;

  cout << "-------------------" << endl;

  DualReducer dr = DualReducer(16, &A, &bl, &bu, &c, &l, &u);
  cout << solMessage(dr.status) << " " << dr.exe_solve << endl;
  double r0_obj = dr.duals[0]->score;
  double gap = (r0_obj - dr.best_score) / r0_obj * 100;
  cout << dr.best_score << " " << r0_obj << " " << gap << "%" << endl;
  
  // testL3Cache();

  // int n = 2;
  // int m = 4;
  // MatrixXd A (m, n); A << 1.0/3, 1, 4.0/3, -1, -1.0/4, 1, 1, -1;
  // VectorXd bl (m); bl << 8.0/3, -DBL_MAX, 1.0/2, 0;
  // VectorXd bu (m); bu << DBL_MAX, 17.0/3, 3, DBL_MAX;
  // VectorXd c (n); c << 1, 1;
  // VectorXd l (n); l << 0, 0;
  // VectorXd u (n); u << 10, 10;

  // Dual dual = Dual(16, A, bl, bu, c, l, u);
  // cout << solMessage(dual.status) << endl;
  // cout << dual.exe_solve << " " << dual.iteration_count << " " << dual.mini_iteration_count << endl;
  // cout << dual.score << endl;

  // GurobiSolver gs = GurobiSolver(A, bl, bu, c, l, u);
  // gs.solveRelaxed();
  // cout << gs.exe_relaxed << " " << gs.iteration_count << endl;
  // cout << gs.relaxed_cscore << endl;
  // // print(dual.x);
  // // print(gs.r0);
  // if (gs.relaxed_status == LS_FOUND && !isEqual(gs.relaxed_cscore, dual.score, 1e-8)){
  //   for (int i = 0; i < n; i ++){
  //     assert(isEqual(dual.x(i), gs.r0(i), 1e-6));
  //   }
  // }
}

int main(){
  //quickRun();
  centroid_test();
}