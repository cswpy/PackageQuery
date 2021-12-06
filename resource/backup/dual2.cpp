#define FMT_HEADER_ONLY

#include <float.h>
#include <iostream>
#include <boost/dynamic_bitset.hpp>
#include "utility.h"
#include "dual.h"
#include "parallel_pq.h"
#include "fmt/core.h"

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

  // Obj
  // Z = 0;

  // Primal solution
  x.resize(n+m); x.fill(0);

  // Leaving variable index in the basis
  int r;
  // Leaving variable index
  int p;
  // Entering variable index
  int q;

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
  vector<pair<double, int>>* init;
  // init size
  int init_size = 0;
  // Parallel PQ
  ParallelPQ* ppq;
  // Dual step
  double dual_step = 0;
  // Primal step
  double primal_step = 0;

  // alpha_q
  VectorXd alpha_q (m);
  // tau for DSE
  VectorXd tau (m);
  
  // Tau
  vector<int> Tau;
  // Delta Z
  // double delta_Z = 0;
  // Tilde a
  VectorXd tilde_a (m);

  cout << "A" << endl;
  cout << A << endl;
  cout << "d" << endl;
  print(d);
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
    // #pragma omp barrier
    // double local_Z = 0;
    // #pragma omp for nowait
    // for (int i = 0; i < n; i ++){
    //   local_Z += x(i)*c(i);
    // }
    // #pragma omp atomic
    // Z += local_Z;
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
        cout << "X" << endl;
        print(x);
        cout << "BOUNDS" << endl;
        print(l);
        print(u);
        print(bl);
        print(bu);
        cout << "END" << endl;
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

        double local_init_size = 0;
        double local_init_start = 0;
        vector<int> local_init ((n+m)/core+1);

        #pragma omp barrier
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
          init = new vector<pair<double, int>>(init_size);
        }
        #pragma omp barrier
        for (int i = local_init_start; i < local_init_start + local_init_size; i ++){
          int j = local_init[i-local_init_start];
          (*init)[i] = {-d(j)/alpha_r(j), j};
        }
        pro.stop(8);
      } else{
        double local_score = 0;
        #pragma omp for nowait
        for (int i = 0; i < n; i ++) local_score += x(i) * cc(i);
        #pragma omp atomic
        score += local_score;
      }
    }
    if (status != LS_FOUND){
      // Ratio test using Algorithm 4
      pro.clock(9, false);
      ppq = new ParallelPQ(core, init);
      pro.stop(9, false);
      q = -1;
      pro.clock(10, false);
      int skip_count = 0;
      cout << "ALPHA_R" << endl;
      print(alpha_r);
      while (ppq->size() > 0){
        q = ppq->peak().second;

        
        double reduction;
        if (q < n) reduction = (u(q)-l(q))*fabs(alpha_r(q))*sign_max_delta;
        else{
          int index = q - n;
          reduction = (bu(index)-bl(index))*fabs(alpha_r(q))*sign_max_delta;
        }
        cout << iteration_count << " " << ppq->size() << " q:" << q << " " << max_delta << " " << alpha_r(q) << " RED " << reduction << endl;

        if (isGreaterEqual(fabs(max_delta), fabs(reduction))){
          max_delta -= reduction;
          ppq->pop();
          skip_count ++;
        } else break;
      }
      if (ppq->size() == 0){
        status = LS_DUAL_UNBOUNDED;
        delete ppq;
        break;
      } else dual_step = d(q) / alpha_r(q) * sign_max_delta;
      cout << "SKIP " << ppq->size() << " "<< skip_count << endl;
      delete ppq;
      pro.stop(10, false);
      pro.clock(11, false);

      // Iteration count for pivot
      iteration_count ++;

      // FTran
      if (q < n) alpha_q = Binv * A.col(q);
      else alpha_q = Binv.col(q-n);

      // DSE FTran
      tau = Binv * rho_r;

      // Basis change and update
      tilde_a.fill(0);
      // delta_Z = 0;
      Tau.clear();
      pro.stop(11, false);
      #pragma omp parallel num_threads(core)
      {
        pro.clock(12);
        // update d using Algorithm 5
        VectorXd local_tilde_a (m); local_tilde_a.fill(0);
        // double local_delta_Z = 0;
        #pragma omp for nowait
        for (int i = 0; i < n+m; i ++){
          if (!inv_bhead[i]){
            d(i) -= dual_step * alpha_r(i) * sign_max_delta;
            if (i < n){
              if (!isEqual(l(i), u(i))){
                if (isEqual(x(i), l(i)) && isLess(d(i), 0)){
                  #pragma omp critical
                  {
                    Tau.push_back(i);
                  }
                  local_tilde_a += (u(i)-l(i)) * A.col(i);
                  // local_delta_Z += (u(i)-l(i)) * c(i);
                } else if (isEqual(x(i), u(i)) && isGreater(d(i), 0)){
                  #pragma omp critical
                  {
                    Tau.push_back(i);
                  }
                  local_tilde_a -= (u(i)-l(i)) * A.col(i);
                  // local_delta_Z -= (u(i)-l(i)) * c(i);
                }
              }
            } else{
              int index = i-n;
              if (!isEqual(bl(index), bu(index))){
                if (isEqual(x(i), bl(index)) && isLess(d(i), 0)){
                  #pragma omp critical
                  {
                    Tau.push_back(i);
                  }
                  local_tilde_a(index) += (bu(index)-bl(index));
                } else if (isEqual(x(i), bu(index)) && isGreater(d(i), 0)){
                  local_tilde_a(index) -= (bu(index)-bl(index));
                  #pragma omp critical
                  {
                    Tau.push_back(i);
                  }
                }
              }
            }
          }
        }
        for (int i = 0; i < m; i ++){
          #pragma omp atomic
          tilde_a(i) += local_tilde_a(i);
        }
        // #pragma omp atomic
        // delta_Z += local_delta_Z;
        pro.stop(12);
        #pragma omp barrier
        pro.clock(13);
        #pragma omp master
        {
          d(p) = -dual_step;
          if (Tau.size() > 0){
            cout << "Binv" << endl;
            cout << Binv << endl;
            cout << "Tilde_a" << endl;
            print(tilde_a);
            VectorXd delta_xB = Binv * tilde_a;
            cout << "DELTA XB" << endl;
            print(delta_xB);
            for (int i = 0; i < m; i ++){
              x(bhead(i)) -= delta_xB(i);
              // if (bhead(i) < n) delta_Z -= c(bhead(i)) * delta_xB(i);
            }
          }

          cout << "AFT TAU" << endl;
          print(x);

          // Z += delta_Z;

          // Update primal solution
          primal_step = max_delta / alpha_q(r);
          for (int i = 0; i < m; i ++) x(bhead(i)) -= primal_step * alpha_q(i);
          x(q) += primal_step;

          cout << "MAX_DELTA: " << max_delta << " PRIMAL_STEP: " << primal_step << endl;
          print(alpha_q);
          cout << "AFT PRIMAL" << endl;
          print(x);

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

          // Update Z
          // Z += dual_step * delta;
        }
        pro.stop(13);
        // Flip bounds using Algorithm 6
        pro.clock(14);
        #pragma omp for nowait
        for (int i = 0; i < Tau.size(); i ++){
          int j = Tau[i];
          if (j < n){
            if (isEqual(x(j), l(j))) x(j) = u(j);
            else x(j) = l(j);
          } else{
            int index = j - n;
            if (isEqual(x(j), bl(index))) x(j) = bu(index);
            else x(j) = bl(index);
          }
        }
        pro.stop(14);
      }
    } else break;
  }
  //pro.print();
  exe_solve = chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - start).count() / 1000000.0;
}