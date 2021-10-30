#define FMT_HEADER_ONLY

#include <float.h>
#include <iostream>
#include <boost/dynamic_bitset.hpp>
#include "utility.h"
#include "dual.h"
#include "parallel_pq.h"
#include "fmt/core.h"
#include "map_sort.h"

using namespace std;

// Primal feasibility tolerance
const double kE_P = 1e-7;
// Relative primal feasibility tolerance
const double kE_r = 1e-9;
// Dual feasibility tolerance
const double kE_D = 1e-7;
// Pivot tolerance
const double kE_ap = 1e-6;
// Zero tolerance
const double kE_z = 1e-12;
// Drop tolerance
const double kE_0 = 1e-14;

double primalInfeasibilities(double x, double l, double u){
  if (x < l) return x-l;
  if (x > u) return x-u;
  return 0;
}

double getSlope(const VectorXd& u, const VectorXd& l, const VectorXd& bu, const VectorXd& bl, const VectorXd& alpha_r, int i){
  int n = u.size();
  if (i < n) return (u(i)-l(i)) * fabs(alpha_r(i));
  else{
    int index = i-n;
    return (bu(index)-bl(index)) * fabs(alpha_r(i));
  }
}

// All slopes are non-negative
// We want to keep k smallest scores such that the sum is at least max_slope
void pushHeap(const VectorXd& u, const VectorXd& l, const VectorXd& bu, const VectorXd& bl, const VectorXd& alpha_r, const VectorXd& d, pair<double, int>* heap, int &sz, double& sum_slope, double max_slope, int i){
  double score = d(i) / alpha_r(i);
  double current_slope = getSlope(u, l, bu, bl, alpha_r, i);
  if (sum_slope <= max_slope){
    heap[sz] = {score, i};
    sz ++;
    push_heap(heap, heap+sz);
    sum_slope += current_slope;
  } 
  else{
    // For sure n > 0 since sum_slope > max_slope >= 0
    if (score <= heap[0].first){
      double replaced_slope = getSlope(u, l, bu, bl, alpha_r, heap[0].second);
      double next_sum_slope = sum_slope - replaced_slope + current_slope;
      if (next_sum_slope < max_slope){
        heap[sz] = {score, i};
        sz ++;
        push_heap(heap, heap+sz);
        sum_slope += current_slope;
      } else{
        pop_heap(heap, heap+sz);
        heap[sz] = {score, i};
        push_heap(heap, heap+sz);
        sum_slope = next_sum_slope;
      }
    }
  }
}

Dual::~Dual(){
}

// Maximizing cc 
// Dual minimizing c
Dual::Dual(int core, const MatrixXd& AA, const VectorXd& bbl, const VectorXd& bbu,  const VectorXd& cc, const VectorXd& ll, const VectorXd& uu){
  vector<string> names = {"Copy", "RestInit", "BoundStricten", "InitXBound", "Compute XB", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"};
  Profiler pro = Profiler(names);
  chrono::high_resolution_clock::time_point start;
  start = chrono::high_resolution_clock::now();
  status = LS_NOT_FOUND;
  iteration_count = 0;
  score = 0;

  n = cc.size();
  m = bbl.size();

  pro.clock(0);
  MatrixXd A = -AA;
  VectorXd bl = bbl;
  VectorXd bu = bbu;
  VectorXd c = -cc;
  VectorXd u = uu;
  VectorXd l = ll;
  pro.stop(0);
  pro.clock(1);
  // Norm square steepest edges
  VectorXd beta (m); beta.fill(1);

  // Reduced cost
  VectorXd d (n+m); d.fill(0); 
  d(seqN(0, n)) = c;

  // Basis indices
  bhead.resize(m);
  for (int i = 0; i < m; i ++) bhead[i] = i+n;
  // Inv Basis indices
  boost::dynamic_bitset<> inv_bhead {n+m, 0};
  for (int i = 0; i < m; i ++) inv_bhead[bhead[i]].flip();

  // Primal solution
  x.resize(n+m); x.fill(0);

  // Leaving variable index in the basis
  int r;
  // Leaving variable index
  int p;
  // Entering variable index
  int q, q_index;

  // Binv
  MatrixXd Binv = MatrixXd::Identity(m, m);
  // rho_r
  VectorXd rho_r (m);
  // alpha_r
  VectorXd alpha_r (n+m);
  // delta of alpha_r
  double max_delta = 0;
  int sign_max_delta = 0;
  
  // init for BRFT PPQ
  pair<double, int>* init;
  // init size
  int init_size = 0;
  // Dual step
  double dual_step = 0;
  // Primal step
  double primal_step = 0;

  // alpha_q
  VectorXd alpha_q (m);
  // tau for DSE
  VectorXd tau (m);
  
  // slope_partition
  VectorXi slope_partition (core+1);
  // slope_sums_partition
  VectorXd slope_sums_parition (core);
  // slope_sums
  VectorXd slope_sums;

  // Tilde a
  VectorXd tilde_a (m);

  // cout << "A" << endl;
  // cout << A << endl;
  // cout << "d" << endl;
  // print(d);

  pro.stop(1);
  #pragma omp parallel num_threads(core)
  {
    // Phase-1 Dual
    pro.clock(2);
    // Stricten bound
    #pragma omp master
    {
      for (int i = 0; i < m; i ++){
        if (bbl(i) == -DBL_MAX) bl(i) = 0;
        if (bbu(i) == DBL_MAX) bu(i) = 0;
      }
    }
    #pragma omp barrier
    for (int i = 0; i < m; i ++){
      if (bbl(i) == -DBL_MAX){
        double local_bound = 0;
        #pragma omp for nowait
        for (int j = 0; j < n; j ++){
          if (AA(i, j) > 0) local_bound += ll(j) * AA(i, j);
          else local_bound += uu(j) * AA(i, j);
        }
        #pragma omp atomic
        bl(i) += local_bound;
      }
      if (bbu(i) == DBL_MAX){
        double local_bound = 0;
        #pragma omp for nowait
        for (int j = 0; j < n; j ++){
          if (AA(i, j) > 0) local_bound += uu(j) * AA(i, j);
          else local_bound += ll(j) * AA(i, j);
        }
        #pragma omp atomic
        bu(i) += local_bound;
      }
    }
    pro.stop(2);
    // Initialize x
    pro.clock(3);
    #pragma omp for nowait
    for (int i = 0; i < n; i ++){
      if (d(i) < 0) x(i) = u(i);
      else x(i) = l(i);
    }
    pro.stop(3);
    pro.clock(4);
    #pragma omp barrier
    for (int i = 0; i < m; i ++){
      double local_anxn = 0;
      #pragma omp for nowait
      for (int j = 0; j < n; j ++){
        local_anxn += A(i, j) * x(j);
      }
      #pragma omp atomic
      x(i+n) -= local_anxn;
    }
    pro.stop(4);
  }

  // Phase-2 Dual
  while (true){
    #pragma omp parallel num_threads(core)
    {
      // Pricing
      pro.clock(5);
      #pragma omp master
      {
        r = -1;
        double max_dse = -DBL_MAX;
        double delta = 0;

        cout << "BASIS" << endl;
        print(bhead);
        // cout << "X" << endl;
        // print(x);
        // cout << "BOUNDS" << endl;
        // print(l);
        // print(u);
        // cout << "SLACK BOUNDS" << endl;
        // print(bl);
        // print(bu);
        // cout << "END" << endl;

        for (int i = 0; i < m; i ++){
          if (bhead(i) < n){
            // n structural variables
            int index = bhead(i);
            delta = primalInfeasibilities(x(bhead(i)), l(index), u(index));
          } else{
            // m slack variables
            int index = bhead(i) - n;
            delta = primalInfeasibilities(x(bhead(i)), bl(index), bu(index));
          }
          if (!isEqual(delta, 0)){
            double dse = delta*delta/beta(i);
            if (max_dse < dse){
              max_dse = dse;
              max_delta = delta;
              r = i;
            }
          }
        }
        if (r == -1){
          // Optimal basis found
          status = LS_FOUND;
        } else{
          sign_max_delta = sign(max_delta);
          p = bhead(r);
          // BTran
          rho_r = Binv.row(r);
        }
      }
      pro.stop(5);
      #pragma omp barrier
      if (status != LS_FOUND){
        pro.clock(6);
        // Pivot row
        #pragma omp for
        for (int i = 0; i < n+m; i ++){
          double val = 0;
          if (i < n){
            for (int j = 0; j < m; j ++) 
              val += A(j, i) * rho_r(j);
          } else{
            val += rho_r(i-n);
          }
          alpha_r(i) = val * sign_max_delta;
        }
        pro.stop(6);
        pro.clock(7);
        // Initialize Ratio test
        #pragma omp master
        {
          init_size = 0;
        }
        #pragma omp barrier
        int local_init_size = 0;
        int local_init_start = 0;
        if (iteration_count == 0){
          vector<int> local_init ((n+m)/core+1);
          #pragma omp for nowait
          for (int i = 0; i < n+m; i ++){
            if (!inv_bhead[i]){
              // Non-basic
              if (i < n){
                if ((isGreater(alpha_r(i), 0) && isEqual(x(i), l(i))) || (isLess(alpha_r(i), 0) && isEqual(x(i), u(i)))){
                  local_init[local_init_size] = i;
                  local_init_size ++;
                }
              } else{
                if ((isGreater(alpha_r(i), 0) && isEqual(x(i), bl(i-n))) || (isLess(alpha_r(i), 0) && isEqual(x(i), bu(i-n)))){
                  local_init[local_init_size] = i;
                  local_init_size ++;
                }
              }
            }
          }
          #pragma omp critical
          {
            local_init_start = init_size;
            init_size += local_init_size;
          }
          pro.stop(7);
          pro.clock(8);
          #pragma omp barrier
          #pragma omp master
          {
            init = new pair<double, int>[init_size];
          }
          #pragma omp barrier
          for (int i = local_init_start; i < local_init_start + local_init_size; i ++){
            int j = local_init[i-local_init_start];
            init[i] = {d(j)/alpha_r(j), j};
          }
          pro.stop(8);
        } else{
          double sum_slope = 0;
          pair<double, int>* local_init = new pair<double, int>[(n+m)/core+1];
          #pragma omp for nowait
          for (int i = 0; i < n+m; i ++){
            if (!inv_bhead[i]){
              // Non-basic
              if (i < n){
                if ((isGreater(alpha_r(i), 0) && isEqual(x(i), l(i))) || (isLess(alpha_r(i), 0) && isEqual(x(i), u(i)))){
                  pushHeap(u, l, bu, bl, alpha_r, d, local_init, local_init_size, sum_slope, fabs(max_delta), i);
                }
              } else{
                if ((isGreater(alpha_r(i), 0) && isEqual(x(i), bl(i-n))) || (isLess(alpha_r(i), 0) && isEqual(x(i), bu(i-n)))){
                  pushHeap(u, l, bu, bl, alpha_r, d, local_init, local_init_size, sum_slope, fabs(max_delta), i);
                }
              }
            }
          }
          #pragma omp critical
          {
            local_init_start = init_size;
            init_size += local_init_size;
          }
          pro.stop(7);
          pro.clock(8);
          #pragma omp barrier
          #pragma omp master
          {
            init = new pair<double, int>[init_size];
          }
          #pragma omp barrier
          for (int i = local_init_start; i < local_init_start + local_init_size; i ++){
            init[i] = local_init[i - local_init_start];
          }
          delete[] local_init;
          pro.stop(8);
        }
      } else{
        double local_score = 0;
        #pragma omp for nowait
        for (int i = 0; i < n; i ++) local_score += x(i) * cc(i);
        #pragma omp atomic
        score += local_score;
      }
    }
    if (status != LS_FOUND){
      // Ratio test using Algorithm 4 combining Algorithm 5
      // cout << "A " << iteration_count << " " << init_size << endl;
      pro.clock(9, false);
      map_sort::Sort(init, init_size, core);
      pro.stop(9, false);
      slope_partition.fill(init_size);
      slope_sums.resize(init_size);
      #pragma omp parallel num_threads(core)
      {
        pro.clock(10);
        double sum_slope = 0;
        bool is_first = false;
        #pragma omp for
        for (int i = 0; i < init_size; i ++){
          if (!is_first){
            is_first = true;
            slope_partition(omp_get_thread_num()) = i;
          }
          sum_slope += getSlope(u, l, bu, bl, alpha_r, init[i].second);
          slope_sums(i) = sum_slope;
        }
        #pragma omp master
        {
          // cout << "B " << iteration_count << " " << init_size << endl;
          sort(slope_partition.begin(), slope_partition.end());
          slope_sums_parition(0) = 0;
          for (int i = 1; i < core; i ++){
            slope_sums_parition(i) = slope_sums_parition(i-1) + slope_sums(slope_partition(i)-1);
          }
        }
        #pragma omp barrier
        int thread_num = omp_get_thread_num();
        for (int i = slope_partition(thread_num); i < slope_partition(thread_num+1); i ++){
          slope_sums(i) += slope_sums_parition(thread_num);
        }
        #pragma omp barrier
        pro.stop(10);
        pro.clock(11);
        #pragma omp master
        {
          // cout << "ALPHA_R" << endl;
          // print(alpha_r);
          // for (int i = 0; i < init_size; i ++) cout << init[i].first << " " << init[i].second << endl;
          // cout << "MAX_SLOPE: " << fabs(max_delta) << endl;
          // print(slope_sums);
          // cout << "C " << iteration_count << " " << init_size << endl;
          if (isLessEqual(slope_sums(init_size-1), fabs(max_delta))){
            status = LS_DUAL_UNBOUNDED;
          } else{
            q_index = upper_bound(slope_sums.begin(), slope_sums.end(), fabs(max_delta)) - slope_sums.begin();
            q = init[q_index].second;
            if (q_index > 0) max_delta = sign_max_delta * (fabs(max_delta) - slope_sums(q_index-1));
            dual_step = d(q) / alpha_r(q) * sign_max_delta;
            d(p) = -dual_step;
            tilde_a.fill(0);

            // FTran
            if (q < n) alpha_q = Binv * A.col(q);
            else alpha_q = Binv.col(q-n);

            // DSE FTran
            tau = Binv * rho_r;
          }
          // cout << "C1 " << iteration_count << " " << init_size << endl;
        }
        pro.stop(11);
        #pragma omp barrier
        if (status != LS_DUAL_UNBOUNDED){
          pro.clock(12);
          // Update d
          #pragma omp for nowait
          for (int i = 0; i < n+m; i ++){
            if (!inv_bhead[i]){
              // Non-basic
              d(i) -= dual_step * alpha_r(i) * sign_max_delta;
            }
          }
          pro.stop(12);
          pro.clock(13);
          // Update tilde_a
          VectorXd local_tilde_a (m); local_tilde_a.fill(0);
          #pragma omp for nowait
          for (int i = 0; i < q_index; i ++){
            int j = init[i].second;
            if (j < n){
              if (isEqual(x(j), l(j))){
                local_tilde_a += (u(j) - l(j)) * A.col(j);
                x(j) = u(j);
              } else{
                local_tilde_a += (l(j) - u(j)) * A.col(j);
                x(j) = l(j);
              }
            } else{
              int index = j - n;
              if (isEqual(x(j), bl(index))){
                local_tilde_a(index) += (bu(index) - bl(index));
                x(j) = bu(index);
              } else{
                local_tilde_a(index) += (bl(index) - bu(index));
                x(j) = bl(index);
              }
            }
          }
          for (int i = 0; i < m; i ++){
            #pragma omp atomic
            tilde_a(i) += local_tilde_a(i);
          }
          #pragma omp barrier
          pro.stop(13);
          pro.clock(14);
          #pragma omp master
          {
            // cout << "D " << iteration_count << " " << init_size << endl;
            // Update xb based on tilde_a
            VectorXd delta_xB = Binv * tilde_a;
            for (int i = 0; i < m; i ++){
              x(bhead(i)) -= delta_xB(i);
            }

            // Update primal solution
            primal_step = max_delta / alpha_q(r);
            for (int i = 0; i < m; i ++) x(bhead(i)) -= primal_step * alpha_q(i);
            x(q) += primal_step;

            // Update edge norm beta
            for (int i = 0; i < m; i ++){
              if (i != r){
                double ratio = alpha_q(i) / alpha_q(r);
                beta(i) += ratio * (ratio * beta(r) - 2 * tau(i));
              }
            }
            beta(r) /= (alpha_q(r) * alpha_q(r));

            // Update bhead and inv_bhead
            inv_bhead[p].flip();
            bhead(r) = q;
            inv_bhead[q].flip();

            // Update Binv
            for (int i = 0; i < m; i ++){
              if (i != r){
                double ratio = alpha_q(i) / alpha_q(r);
                Binv.row(i) -= ratio * Binv.row(r);
              }
            }
            Binv.row(r) /= alpha_q(r);

            // Iteration count for pivot
            iteration_count ++;
          }
          pro.stop(14);
        }
      }
      if (status == LS_DUAL_UNBOUNDED) break;
    } else break;
  }
  pro.print();
  exe_solve = chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - start).count() / 1000000.0;
}