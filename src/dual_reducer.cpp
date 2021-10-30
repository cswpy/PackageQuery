#define FMT_HEADER_ONLY

#include <iostream>
#include <float.h>
#include <chrono>
#include <random>
#include <unordered_map>
#include <numeric>

#include "fmt/core.h"
#include "pseudo_walker.h"
#include "gurobi_solver.h"
#include "simplex.h"
#include "utility.h"
#include "dual_reducer.h"
#include "dual.h"
#include "omp.h"

using namespace Eigen;
using namespace std;

const int kIlpThreshold = 50;
const double kStepTolerance = 1.05;

DualReducer::~DualReducer(){
  if (layer_count > 0) delete original_indices[0];
  for (int i = 1; i < layer_count+1; i ++){
    delete duals[i-1];
    delete As[i];
    delete bls[i];
    delete bus[i];
    delete cs[i];
    delete ls[i];
    delete us[i];
    delete original_indices[i];
  }
}

DualReducer::DualReducer(int core, MatrixXd* AA, VectorXd* bbl, VectorXd* bbu, VectorXd* cc, VectorXd* ll, VectorXd* uu){
  chrono::high_resolution_clock::time_point start;
  start = chrono::high_resolution_clock::now();

  As.push_back(AA);
  bls.push_back(bbl);
  bus.push_back(bbu);
  cs.push_back(cc);
  ls.push_back(ll);
  us.push_back(uu);

  int nn = cc->size();
  int m = bbl->size();

  original_indices.push_back(new VectorXi(nn));
  best_x.resize(nn);
  best_score = 0;
  layer_count = 0;
  status = LS_NOT_FOUND;

  MatrixXd* A;
  VectorXd* bl, *bu, *c, *l, *u;
  VectorXi* original_index;

  while (true){
    A = As[layer_count];
    bl = bls[layer_count];
    bu = bus[layer_count];
    c = cs[layer_count];
    l = ls[layer_count];
    u = us[layer_count];
    original_index = original_indices[layer_count];

    int n = c->size();
    Dual* dual = new Dual(core, *A, *bl, *bu, *c, *l, *u);
    duals.push_back(dual);

    if (n < kIlpThreshold){
      if (layer_count == 0){
        for (int i = 0; i < nn; i ++) (*original_index)(i) = i;
      }
      break;
    }

    if (dual->status == LS_FOUND){
      VectorXd centroid_dir (n); centroid_dir.fill(0);
      VectorXd r0 (n);
      VectorXd stay_r0 (n);
      VectorXd sum_row(m); sum_row.fill(0);
      int stay_count = 0;
      #pragma omp parallel num_threads(core)
      {
        if (layer_count == 0){
          #pragma omp for nowait
          for (int i = 0; i < nn; i ++) (*original_index)(i) = i;
        }
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
      }

      PseudoWalker walker = PseudoWalker(centroid_dir, true, core);
      
      // print(dual->x);
      // print(centroid_dir);
      int step_horizon = (int)ceil(n / log2(n) * kStepTolerance);
      while (true){
        int step = walker.step();
        int index = abs(step)-1;
        if (stay_r0(index) >= 0){
          stay_r0(index) = -1;
          stay_count ++;
        }
        r0(index) += sign(step);

        // cout << "STEP " << step_count << " " << step << " " << index << endl;
        // print(r0);

        if (isGreaterEqual((*l)(index)-r0(index), 1) || isGreaterEqual(r0(index)-(*u)(index), 1) || walker.step_count == step_horizon) break;
      }
      cout << layer_count << " " << walker.step_count << " " << stay_count << endl;

      VectorXi stay_index (stay_count);
      VectorXi* next_original_index = new VectorXi(stay_count);
      original_indices.push_back(next_original_index);
      MatrixXd* nextA = new MatrixXd(m, stay_count);
      As.push_back(nextA);
      VectorXd* nextbl = new VectorXd(*bl);
      bls.push_back(nextbl);
      VectorXd* nextbu = new VectorXd(*bu);
      bus.push_back(nextbu);
      VectorXd* nextc = new VectorXd(stay_count);
      cs.push_back(nextc);
      VectorXd* nextl = new VectorXd(stay_count);
      ls.push_back(nextl);
      VectorXd* nextu = new VectorXd(stay_count);
      us.push_back(nextu);
      int start_stay_count = 0;
      #pragma omp parallel num_threads(core)
      {
        int local_start_index = -1;
        int local_stay_count = 0;
        VectorXi local_stay_index (stay_count);
        VectorXd local_change_b (m); local_change_b.fill(0);
        #pragma omp for nowait
        for (int i = 0; i < n; i ++){
          if (stay_r0(i) < 0){
            local_stay_index(local_stay_count) = i;
            local_stay_count ++;
          } 
          else{
            best_x((*original_index)(i)) = stay_r0(i);
            for (int j = 0; j < m; j ++){
              local_change_b(j) += (*A)(j, i) * stay_r0(i);
            }
          }
        }
        #pragma omp critical
        {
          local_start_index = start_stay_count;
          start_stay_count += local_stay_count;
        }
        for (int i = 0; i < m; i ++){
          if ((*bl)(i) != -DBL_MAX){
            #pragma omp atomic
            (*nextbl)(i) -= local_change_b(i);
          }
          if ((*bu)(i) != DBL_MAX){
            #pragma omp atomic
            (*nextbu)(i) -= local_change_b(i);
          }
        }
        #pragma omp barrier
        for (int i = local_start_index; i < local_start_index + local_stay_count; i ++){
          stay_index(i) = local_stay_index(i-local_start_index);
          (*next_original_index)(i) = (*original_index)(stay_index(i));
        }
        #pragma omp barrier
        #pragma omp for nowait
        for (int i = 0; i < stay_count; i ++){
          for (int j = 0; j < m; j ++){
            (*nextA)(j, i) = (*A)(j, stay_index[i]);
          }
          (*nextc)(i) = (*c)(stay_index[i]);
          (*nextl)(i) = (*l)(stay_index[i]);
          (*nextu)(i) = (*u)(stay_index[i]);
        }
      }
      layer_count ++;
    } else{
      status = dual->status;
      break;
    }
    break;
  }

  if (status == LS_NOT_FOUND){
    GurobiSolver gs = GurobiSolver(*A, *bl, *bu, *c, *l, *u);
    gs.solveIlp();
    if (gs.ilp_status == LS_FOUND){
      for (int i = 0; i < gs.n; i ++) best_x((*original_index)(i)) = gs.x0(i);
      #pragma omp parallel num_threads(core)
      {
        double local_best_score = 0;
        #pragma omp for nowait
        for (int i = 0; i < nn; i ++) local_best_score += best_x(i) * ((*cc)(i));
        #pragma omp atomic
        best_score += local_best_score;
      }
    }
    status = gs.ilp_status;
  }

  exe_solve = chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - start).count() / 1000000.0;
}