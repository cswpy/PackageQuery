#define FMT_HEADER_ONLY

#include <iostream>
#include <float.h>
#include <chrono>
#include <random>
#include <unordered_map>

#include "fmt/core.h"
#include "lattice_solver.h"
#include "pseudo_walker.h"
#include "simplex.h"
#include "utility.h"
#include "omp.h"

using namespace Eigen;
using namespace std;

#define tab(row, col) (simplex->tableau[simplex->numcols*(row)+(col)])

LatticeSolver::~LatticeSolver(){
}

LatticeSolver::LatticeSolver(int core, const MatrixXd& A, const VectorXd& b, const VectorXd& c, const VectorXd& u): A(A), b(b), c(c), u(u){
  chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();
  this->core = core;
  m = b.size();
  n = c.size();
  r0.resize(n);
  fracs.resize(n);
  near_r0.resize(n);
  fscores.resize(m);
  bhead.resize(m);
  for (int i = 0; i < m; i ++) fscores(i) = b(i);
  nbasic.resize(n); fill(nbasic.begin(), nbasic.end(), 0);
  cscore = 0;
  relaxed_cscore = 0;
  best_cscore = -DBL_MAX;
  lattice_dirs = CooSparse(n, m);
  exe_init = chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - start).count() / 1000000.0;
  start = chrono::high_resolution_clock::now();
  Simplex* simplex = new Simplex(core, A, b, c, u);
  exe_relaxed = chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - start).count() / 1000000.0;
  start = chrono::high_resolution_clock::now();
  if (simplex->status == LS_FOUND){
    for (int i = 0; i < m; i ++){
      inv_bhead[simplex->bhead[i]] = i;
      bhead[i] = simplex->bhead[i];
    }
    status = LS_NOT_FOUND;
    #pragma omp parallel num_threads(core)
    {
      double local_cscore = 0;
      double local_relaxed_cscore = 0;
      VectorXd local_fscores (m); local_fscores.fill(0);
      #pragma omp for nowait
      for (int i = 0; i < n; i ++){
        r0(i) = tab(1, i);
        fracs(i) = modf(r0(i), &near_r0(i));
        if (fracs(i) > 0.5) near_r0(i) ++;
        local_cscore += c(i)*near_r0(i);
        local_relaxed_cscore += c(i)*r0(i);
        for (int j = 0; j < m; j ++) local_fscores(j) += A(j, i) *near_r0(i);
      }
      #pragma omp atomic
      cscore += local_cscore;
      #pragma omp atomic
      relaxed_cscore += local_relaxed_cscore;
      for (int i = 0; i < m; i ++){
        #pragma omp atomic
        fscores(i) -= local_fscores(i);
      }
      #pragma omp barrier
      #pragma omp for nowait
      for (int i = 0; i < m+n; i ++){
        if (inv_bhead.find(i) == inv_bhead.end()){
          // Non-basic
          int index = i;
          for (int j = 0; j < m; j ++){
            if (i > bhead[j]) index --;
          }
          double dir = 1.0;
          if (isEqual(tab(1, i), 0)) dir = 1.0;
          if (i < n) nbasic[index] = (i+1)*dir;
          for (int j = 2; j <= m+1; j ++){
            int basic_col = bhead[j-2];
            if (basic_col < n){
              lattice_dirs.addEntry(index, basic_col, -tab(j, i) * dir);
            }
          }
        }
      }
    }
  } else status = simplex->status;
  exe_tableau = chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - start).count() / 1000000.0;
  delete simplex;
}

double LatticeSolver::getObjValue(VectorXd sol){
  return c.dot(sol);
}

bool LatticeSolver::checkFeasibility(VectorXd sol){
  for (int i = 0; i < n; i++) {
    if (isGreater(sol(i), u(i), kFeasTolerance) || isLess(sol(i), 0, kFeasTolerance)) {
      return false;
    }
  }
  for (int i = 0; i < m; i++) {
    if (isGreater(A.row(i).dot(sol), b(i), kFeasTolerance)) {
      return false;
    }
  }
  return true;
}

bool LatticeSolver::checkFeasibilityStep(const VectorXd& current_fscores, const unordered_map<int,int>& bound_map){
  if (bound_map.size() > 0) return false;
  for (int j = 0; j < m; j ++){
    if (isLess(current_fscores(j), 0, kFeasTolerance)) return false;
  }
  return true;
}

bool LatticeSolver::latticeStep(int d, VectorXd& x, VectorXd& current_fscores, double& current_cscore, unordered_map<int,int>& bound_map){
  int i = abs(d) - 1;
  int s = sign(d);
  x(i) += s;
  if (bound_map.find(i) == bound_map.end()){
    int lb_violation = floor(x(i));
    int ub_violation = ceil(x(i) - u(i));
    if (lb_violation < 0) bound_map[i] = lb_violation;
    else if (ub_violation > 0) bound_map[i] = ub_violation;
  } else {
    bound_map[i] += s;
    if (bound_map[i] == 0) bound_map.erase(i);
    if (abs(bound_map[i]) > kBoundTol) return false;
  }
  for (int j = 0; j < m; j ++) current_fscores(j) -= s*A(j, i);
  current_cscore += s*c(i);
  return true;
}

void LatticeSolver::solve(double time_budget){
  chrono::high_resolution_clock::time_point ep = chrono::high_resolution_clock::now();
  try_count = 0;
  avg_step_count = 0;
  #pragma omp parallel num_threads(core)
  {
    // Local declaration
    default_random_engine gen {static_cast<long unsigned int>(time(0))};
    //default_random_engine gen {0};
    uniform_real_distribution u_dist(0.0, 1.0);
    int local_try_count = 0;
    int local_step_count = 0;
    double local_best_cscore = -DBL_MAX;
    VectorXd local_ts (20); local_ts.fill(0);
    VectorXd local_best_x;
    while (true) {
      if (chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - ep).count() / 1000000000.0 > time_budget)
        break;
      local_try_count ++;
      VectorXd coefs (n);
      coefs(0) = 1;
      for (int i = 1; i < n; i ++) coefs(i) = u_dist(gen);
      sort(coefs.begin(), coefs.end());
      for (int i = n-1; i >= 1; i--) coefs(i) -= coefs(i-1);
      VectorXd dir (n); dir.fill(0);
      lattice_dirs.inplaceVectorProduct(coefs, dir);
      for (int i = 0; i < n; i ++){
        int comp = abs(nbasic[i]) - 1;
        if (nbasic[i] != 0) dir(comp) = sign(nbasic[i]) * coefs(i);
      }
      PseudoWalker* walker = new PseudoWalker(dir);
      VectorXd x = near_r0;
      VectorXd current_fscores = fscores;
      double current_cscore = cscore;
      unordered_map<int, int> bound_map;
      for (int i = 0; i < n; i ++){
        if (u_dist(gen) < fracs(i)){
          if (fracs(i) <= 0.5) latticeStep(i+1, x, current_fscores, current_cscore, bound_map);
        } else{
          if (fracs(i) > 0.5) latticeStep(-(i+1), x, current_fscores, current_cscore, bound_map);
        }
      }
      if (checkFeasibilityStep(current_fscores, bound_map)){
        if (local_best_cscore < current_cscore){
          local_best_cscore = current_cscore;
          local_best_x = x;
        }
      }
      int step_count = 0;
      double stubborn_scale = 1.0;
      while (true){
        int d = walker->step();
        if (!latticeStep(d, x, current_fscores, current_cscore, bound_map)) break;
        if (checkFeasibilityStep(current_fscores, bound_map)){
          if (local_best_cscore < current_cscore){
            local_best_cscore = current_cscore;
            local_best_x = x;
            stubborn_scale = 1.0;
          }
          if (u_dist(gen) < stubborn_scale * (best_cscore - current_cscore) / best_cscore) break;
        }
        step_count ++;
        stubborn_scale *= kStubbornMultiplier;
        if (step_count > kSizeTol) break;
      }
      delete walker;
      local_step_count += step_count;
    }
    #pragma omp critical 
    {
      if (best_cscore < local_best_cscore){
        best_cscore = local_best_cscore;
        best_x = local_best_x;
      }
    }
    #pragma omp atomic
    try_count += local_try_count;
    #pragma omp atomic
    avg_step_count += local_step_count;
  }
  avg_step_count /= try_count;
  if (best_cscore != -DBL_MAX){
    if (status == LS_NOT_FOUND) status = LS_FOUND;
  }
  exe_solved = chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - ep).count() / 1000000.0;
}

void LatticeSolver::report(){
  fmt::print("-------------Problem statistics-------------\n");
  fmt::print("Problem size: {}x{}\n", m, n);
  fmt::print("Number of non-zero entries: {}\n", lattice_dirs.getSize());
  fmt::print("Time for initialization: {:.2Lf}ms\n", exe_init);
  fmt::print("Time for finding relaxed solution: {:.2Lf}ms\n", exe_relaxed);
  fmt::print("Time for processing simplex tableau: {:.2Lf}ms\n", exe_tableau);
  fmt::print("Time for lattice walking: {:.2Lf}ms\n", exe_solved);
  fmt::print("Total time: {:.2Lf}ms\n", exe_init+exe_relaxed+exe_tableau+exe_solved);
  fmt::print("Total cores: {}\n", core);
  fmt::print("------------Solution statistics-------------\n");
  fmt::print("Total directions: {}\n", try_count);
  fmt::print("Average number of steps per direction: {:.2Lf}\n", avg_step_count);
  fmt::print("Solution status: {}\n", solMessage(status));
  if (status == LS_FOUND){
    fmt::print("Best found combination: {}\n", solCombination(best_x));
    fmt::print("Best found objective: {:.8Lf}\n", best_cscore);
    fmt::print("Relaxed objective: {:.8Lf}\n", relaxed_cscore);
    double relaxed_gap = (relaxed_cscore - best_cscore) / fabs(relaxed_cscore) * 100;
    fmt::print("Objective gap to relaxed solution: {:.2Lf}%\n", relaxed_gap);
  } else {
    fmt::print("Relaxed objective: {:.8Lf}\n", relaxed_cscore);
  }
}

void LatticeSolver::compareReport(VectorXd sol, double solved_time){
  report();
  if (status == LS_FOUND){
    fmt::print("Comparing combination: {}\n", solCombination(sol));
    double comparing_score = getObjValue(sol);
    fmt::print("Comparing objective: {:.8Lf}\n", comparing_score);
    double comparing_gap = (comparing_score - best_cscore) / fabs(comparing_score) * 100;
    double lp_gap = (relaxed_cscore - comparing_score) / fabs(relaxed_cscore) * 100;
    fmt::print("Objective gap to comparing solution: {:.2Lf}%\n", comparing_gap);
    fmt::print("Comparing gap to relaxed solution: {:.2Lf}%\n", lp_gap);
    if (solved_time > 0) fmt::print("Comparing solution found in: {:.2Lf}ms\n", solved_time);
  }
}