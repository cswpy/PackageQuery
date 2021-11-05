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

const int kIlpThreshold = 500;
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

DualReducer::DualReducer(int core, RMatrixXd* AA, VectorXd* bbl, VectorXd* bbu, VectorXd* cc, VectorXd* ll, VectorXd* uu, int opt){
  // vector<string> names = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10"};
  // Profiler pro = Profiler(names);
  chrono::high_resolution_clock::time_point start;
  start = chrono::high_resolution_clock::now();
  // pro.clock(0, false);
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

  RMatrixXd* A;
  VectorXd* bl, *bu, *c, *l, *u;
  VectorXi* original_index;
  // pro.stop(0, false);
  while (true){
    // cout << "Z " << layer_count << endl;
    A = As[layer_count];
    bl = bls[layer_count];
    bu = bus[layer_count];
    c = cs[layer_count];
    l = ls[layer_count];
    u = us[layer_count];
    original_index = original_indices[layer_count];
    // cout << "Z1 " << layer_count << endl;
    int n = c->size();
    if (n < kIlpThreshold){
      if (layer_count == 0){
        // cout << "Z11 " << layer_count << endl;
        for (int i = 0; i < nn; i ++) (*original_index)(i) = i;
        // pro.clock(1, false);
        Dual* dual = new Dual(core, *A, *bl, *bu, *c, *l, *u);
        // pro.stop(1, false);
        duals.push_back(dual);
        // cout << "Z12 " << layer_count << endl;
      }
      break;
    }
    // cout << "Z2 " << layer_count << endl;
    // pro.clock(1, false);
    Dual* dual = new Dual(core, *A, *bl, *bu, *c, *l, *u);
    // pro.stop(1, false);
    duals.push_back(dual);
    //cout << "A " << layer_count << endl;

    if (dual->status == LS_FOUND){
      VectorXd norms (n+m); norms.fill(0);
      VectorXd centroid_dir (n); centroid_dir.fill(0);
      VectorXd r0 (n);
      VectorXd stay_r0 (n);
      VectorXd sum_row(m); sum_row.fill(0);
      RMatrixXd BinvA (m, n+m);
      int stay_count = 0;
      bool is_conservative = true;
      int horizon = (int)ceil(n / log2(n) * kStepTolerance);
      #pragma omp parallel num_threads(core)
      {
        // pro.clock(2);
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
          if (isEqual(frac, 0) && whole == 0) stay_r0(i) = whole;
          else{
            stay_r0(i) = -1;
            local_stay_count ++;
          }
        }
        #pragma omp atomic
        stay_count += local_stay_count;
        // #pragma omp barrier
        // #pragma omp master
        // {
        //   if (stay_count > horizon) is_conservative = false;
        //   stay_count = 0;
        // }
        // #pragma omp barrier
        // if (!is_conservative){
        //   #pragma omp for nowait
        //   for (int i = 0; i < n; i ++){
        //     if (stay_r0(i) == -1){
        //       frac = modf(dual->x(i), &whole);
        //       if (isEqual(frac, 0)){
        //         stay_r0(i) = whole;
        //         local_stay_count --;
        //       }
        //     }
        //   }
        //   #pragma omp atomic
        //   stay_count += local_stay_count;
        // }
        // for (int j = 0; j < m; j ++){
        //   #pragma omp for nowait
        //   for (int i = 0; i < n+m; i ++){
        //     if (!dual->inv_bhead[i]){
        //       if (i < n){
        //         double val = 0;
        //         for (int k = 0; k < m; k ++) val += dual->Binv(j, k) * dual->A(k, i);
        //         BinvA(j, i) = val;
        //       } else{
        //         int index = i - n;
        //         BinvA(j, i) = dual->Binv(j, index);
        //       }
        //     }
        //   }
        // }
        // #pragma omp barrier
        // #pragma omp master
        // {
        //   for (int i = 0; i < m; i ++){
        //     int basic_col = dual->bhead(i);
        //     if (basic_col < n){
        //       //cout << basic_col << " " << dual->x(basic_col) << " " << dual->l(basic_col) << " " << dual->u(basic_col) << endl;
        //       if (isEqual(dual->x(basic_col), dual->l(basic_col)) || isEqual(dual->x(basic_col), dual->u(basic_col))){
        //         cout << "DENGERACY DETECTED AT " << basic_col << endl;
        //       }
        //     } else{
        //       int index = basic_col - n;
        //       //cout << basic_col << " " << dual->x(basic_col) << " " << dual->bl(index) << " " << dual->bu(index) << endl;
        //       if (isEqual(dual->x(basic_col), dual->bl(index)) || isEqual(dual->x(basic_col), dual->bu(index))){
        //         cout << "DENGERACY DETECTED AT " << basic_col << endl;
        //       }
        //     }
        //   }
        // }
        // #pragma omp barrier
        if (opt == 0){
          #pragma omp for nowait
          for (int i = 0; i < m+n; i ++){
            if (!dual->inv_bhead[i]){
              double norm = 0;
              if (i < n){
                norm ++;
                for (int j = 0; j < m; j ++){
                  if (dual->bhead(j) < n){
                    double val = 0;
                    for (int k = 0; k < m; k ++){
                      val += dual->Binv(j, k) * dual->A(k, i);
                    }
                    norm += val*val;
                  }
                }
              } else{
                int index = i - n;
                for (int j = 0; j < m; j ++){
                  if (dual->bhead(j) < n){
                    norm += dual->Binv(j, index)*dual->Binv(j, index);
                  }
                }
              }
              norms(i) = sqrt(norm);
            }
          }
          #pragma omp barrier
          for (int j = 0; j < m; j ++){
            double local_sum_row = 0;
            #pragma omp for nowait
            for (int i = 0; i < m+n; i ++){
              if (!dual->inv_bhead[i]){
                // Non-basic
                if (i < n){
                  if (isEqual(dual->x(i), dual->l(i))) local_sum_row += dual->A(j, i) / norms(i);
                  else local_sum_row -= dual->A(j, i) / norms(i);
                } else{
                  int index = i - n;
                  if (index == j){
                    if (isEqual(dual->x(i), dual->bl(index))) local_sum_row -= 1.0 / norms(i);
                    else local_sum_row += 1.0 / norms(i);
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
              if (isEqual(dual->x(i), dual->l(i))) centroid_dir(i) = 1.0/norms(i)/n;
              else centroid_dir(i) = -1.0/norms(i)/n;
            }
          }
          #pragma omp barrier
          #pragma omp master
          {
            //cout << "A1 " << layer_count << endl;
            for (int i = 0; i < m; i ++){
              if (dual->bhead(i) < n){
                double val = 0;
                for (int j = 0; j < m; j ++){
                  val += dual->Binv(i, j) * sum_row(j);
                }
                //cout << "HERE " << i << " " << dual->bhead(i) << " VAL:" << val << endl;
                centroid_dir(dual->bhead(i)) = val / n;
              }
            }
          }
          // pro.stop(2);
        }
      }
      // pro.clock(3, false);
      //cout << "A2 " << layer_count << endl;
      //cout << is_conservative << " " << stay_count << endl;
      if (opt == 0){
        PseudoWalker walker = PseudoWalker(centroid_dir, true, core);
        // pro.stop(3, false);

        // print(dual->x);
        // print(centroid_dir);
        // RMatrixXd BA = dual->Binv * dual->A;
        // cout << BA << endl;
        // cout << endl;
        // cout << dual->Binv << endl;
        // cout << "-----------------------" << endl;
        // print(r0);
        // print(centroid_dir);
        
        // pro.clock(4, false);
        while (true){
          int step = walker.step();
          int index = abs(step)-1;
          if (stay_r0(index) >= 0){
            stay_r0(index) = -1;
            stay_count ++;
          }
          r0(index) += sign(step);

          // cout << "STEP " << walker.step_count << " " << sign(step) << " " << index << endl;
          // shortPrint(r0);
          // for (int i = 0; i < n; i ++){
          //   if (!isEqual(r0(i),0)) cout << i << ":" << centroid_dir(i) << " ";
          // }
          // cout << endl;

          // if (stay_count == n/2){
          //   break;
          // }

          if (walker.step_count == horizon){
            break;
          }

          // if (((*l)(index)-r0(index) >= 1) || (r0(index)-(*u)(index) >= 1) || walker.step_count == horizon){
          //   break;
          // }
        }
        // cout << "LAYER: " << layer_count << " STEP COUNT: " << walker.step_count << " STAY:" << stay_count << endl;
      }

      // pro.stop(4, false);
      // pro.clock(5, false);
      VectorXi stay_index (stay_count);
      VectorXi* next_original_index = new VectorXi(stay_count);
      original_indices.push_back(next_original_index);
      RMatrixXd* nextA = new RMatrixXd(m, stay_count);
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
      // pro.stop(5, false);
      //cout << "B " << layer_count << endl;
      #pragma omp parallel num_threads(core)
      {
        // pro.clock(6);
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
        // #pragma omp single
        //cout << "B1 " << layer_count << endl;
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
        // pro.stop(6);
        // pro.clock(7);
        #pragma omp barrier
        for (int i = local_start_index; i < local_start_index + local_stay_count; i ++){
          stay_index(i) = local_stay_index(i-local_start_index);
          (*next_original_index)(i) = (*original_index)(stay_index(i));
        }
        #pragma omp barrier
        #pragma omp for nowait
        for (int i = 0; i < stay_count; i ++){
          for (int j = 0; j < m; j ++){
            (*nextA)(j, i) = (*A)(j, stay_index(i));
          }
          (*nextc)(i) = (*c)(stay_index(i));
          (*nextl)(i) = (*l)(stay_index(i));
          (*nextu)(i) = (*u)(stay_index(i));
        }
        // pro.stop(7);
        // #pragma omp single
        //cout << "B2 " << layer_count << endl;
      }
      layer_count ++;
      if (stay_count >= horizon){
        A = As[layer_count];
        bl = bls[layer_count];
        bu = bus[layer_count];
        c = cs[layer_count];
        l = ls[layer_count];
        u = us[layer_count];
        original_index = original_indices[layer_count];
        break;
      }
      //cout << "B3 " << layer_count << endl;
    } else{
      status = dual->status;
      break;
    }
  }
  //cout << "C " << layer_count << endl;
  //cout << bl->size() << " " << A->innerSize() << " " << A->outerSize() << endl;
  // pro.clock(8, false);
  if (c->size() == 0){
    status = LS_FOUND;
    double local_best_score = 0;
    #pragma omp for nowait
    for (int i = 0; i < nn; i ++) local_best_score += best_x(i) * ((*cc)(i));
    #pragma omp atomic
    best_score += local_best_score;
  }
  if (status == LS_NOT_FOUND){
    //cout << "C01 " << layer_count << endl;
    GurobiSolver gs = GurobiSolver(*A, *bl, *bu, *c, *l, *u);
    //cout << "C02 " << layer_count << endl;
    gs.solveIlp();
    //cout << "C03" << layer_count << endl;
    if (gs.ilp_status == LS_FOUND){
      for (int i = 0; i < gs.n; i ++) best_x((*original_index)(i)) = gs.x0(i);
      //cout << "C1 " << layer_count << endl;
      #pragma omp parallel num_threads(core)
      {
        double local_best_score = 0;
        #pragma omp for nowait
        for (int i = 0; i < nn; i ++) local_best_score += best_x(i) * ((*cc)(i));
        #pragma omp atomic
        best_score += local_best_score;
      }
      //cout << "C2 " << layer_count << endl;
    }
    status = gs.ilp_status;
  }
  // pro.stop(8, false);
  // pro.print();
  exe_solve = chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - start).count() / 1000000.0;
}