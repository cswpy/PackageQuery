#define FMT_HEADER_ONLY

#include <iostream>
#include <float.h>
#include <chrono>
#include <random>
#include <unordered_map>

#include "fmt/core.h"
#include "pseudo_walker.h"
#include "simplex.h"
#include "utility.h"
#include "reducer.h"
#include "omp.h"

using namespace Eigen;
using namespace std;

#define tab(row, col) (simplex.tableau[simplex.numcols*(row)+(col)])

Reducer::~Reducer(){
}

Reducer::Reducer(int core, MatrixXd* AA, VectorXd* bb, VectorXd* cc, VectorXd* uu): A(AA), b(bb), c(cc), u(uu){
  while (true){
    int m = b->size();
    int n = c->size();
    if (n < kReduce) break;
    Simplex simplex = Simplex(core, *A, *b, *c, *u);
    if (simplex.status == LS_FOUND){
      VectorXd centroid_dir (n);
      VectorXd r0 (n);
      unordered_map<int, int> inv_bhead;
      for (int i = 0; i < m; i ++) inv_bhead[simplex.bhead[i]] = i;
      #pragma omp parallel num_threads(core)
      {
        double whole;
        #pragma omp for
        for (int i = 0; i < n; i ++){
          centroid_dir(i) = 0;
          whole = round(tab(1, i));
          if (isEqual(whole-tab(1, i), 0)) r0(i) = whole;
          else r0(i) = -1;
        }
        unordered_map<int, double> local_centroid_dir;
        #pragma omp for nowait
        for (int i = 0; i < m+n; i ++){
          if (inv_bhead.find(i) == inv_bhead.end()){
            // Non-basic
            int index = i;
            for (int j = 0; j < m; j ++){
              if (i > simplex.bhead[j]) index --;
            }
            double dir = -1.0;
            if (isEqual(tab(1, i), 0)) dir = 1.0;
            if (i < n) local_centroid_dir[i] += (dir/n);
            for (int j = 2; j <= m+1; j ++){
              int basic_col = simplex.bhead[j-2];
              if (basic_col < n) local_centroid_dir[basic_col] -= (tab(j, i) * dir / n);
            }
          }
        }
        for (const auto& pi : local_centroid_dir){
          #pragma omp atomic
          centroid_dir(pi.first) += pi.second;
        }
      }
      PseudoWalker walker = PseudoWalker(centroid_dir, true, core);
      int step_horizon = (int)ceil(n / log2(n) * kRedundant);
      for (int i = 0; i < step_horizon; i ++){
        int index = abs(walker.step())-1;
        r0(index) = -1;
      }
      showHistogram(r0, 10, 0, 0);
    } else{
      status = simplex.status;
      break;
    }
    break;
  }
}