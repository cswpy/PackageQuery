#include "pb/util/udeclare.h"
#include "pb/util/uconfig.h"
#include "pb/util/udebug.h"
#include "pb/det/det_prob.h"
#include "pb/core/dual_reducer.h"
#include "pb/core/gurobi_solver.h"
#include "pb/det/det_bound.h"
#include "pb/det/dlv.h"
#include "pb/det/synthetic.h"
#include "pb/det/lsr_prob.h"
#include "pb/util/upostgres.h"
#include "pb/det/lsr.h"
#include "pb/det/det_sql.h"
#include "pb/core/dist.h"
#include "pb/det/det_exp.h"

#include "pb/lib/praxis.hpp"
#include "pb/lib/normal.hpp"
#include "pb/lib/log_normal_truncated_ab.hpp"
#include "pb/lib/log_normal.hpp"
#include "pb/lib/truncated_normal.hpp"
#include "pb/lib/toms178.hpp"
#include "pb/lib/brent.hpp"
#include "pb/core/vg.h"

#include "pb/lib/random_quantile.h"
#include "pb/lib/common.h"

using namespace pb;



void dlv_binary_test(){
  // double left = 0;
  // double dif = samples[N-1] - samples[0];
  // double right = 0.25*dif*dif; // Popoviciu inequality
  // while (abs(left-right) > 1e-8){
  //   double mid = (left+right)/2;
  //   int cnt = 1;
  //   ScalarMeanVar smv = ScalarMeanVar();
  //   for (int i = 0; i < N; i ++){
  //     smv.add(samples[i]);
  //     if (smv.getBiasedVar() > mid){
  //       smv.reset();
  //       smv.add(samples[i]);
  //       cnt ++;
  //     }
  //   }
  //   // cout << cnt << " " << g << " " << left << " " << mid << " " << right << endl;
  //   if (cnt < g) right = mid;
  //   else left = mid;
  // }
  // ScalarMeanVar smv_opt_var = ScalarMeanVar();
  // double opt_var = (left+right) / 2;
  // vector<int> opt_sizes (g, 0);
  // int index = 0;
  // ScalarMeanVar smv = ScalarMeanVar();
  // for (int i = 0; i < N; i ++){
  //   ScalarMeanVar tmp_smv = smv;
  //   smv.add(samples[i]);
  //   if (smv.getBiasedVar() > opt_var){
  //     smv_opt_var.add(tmp_smv.getBiasedVar());
  //     smv.reset();
  //     smv.add(samples[i]);
  //     index ++;
  //   }
  //   opt_sizes[index] ++;
  // }
  // smv_opt_var.add(smv.getBiasedVar());
}

void dlv_graph_test(){
  DetExp exp = DetExp("DLV2");
  int N = 100000;
  VectorXd samples;
  double mean = 0;
  double var = 1;
  double a = mean - sqrt(3*var);
  double b = mean + sqrt(3*var);
  UNUSED(a);
  UNUSED(b);
  Normal dist = Normal(mean, var);
  // Uniform dist = Uniform(a, b);
  dist.sample(samples, N);
  sort(samples.begin(), samples.end());
  ScalarMeanVar whole;
  for (int i = 0; i < (int) samples.size(); i ++){
    whole.add(samples[i]);
  }
  double sigmas = whole.getBiasedVar();
  double upper = 0.25*(whole.getMax() - whole.getMin())*(whole.getMax() - whole.getMin());
  double lower = 0;
  int cur_count = N;

  while (cur_count > 1){
    double left = lower;
    double right = upper;
    double my_cnt, my_vRAC, my_vRWC;
    while (abs(right-left)>1e-6){
      double mid = (left+right)/2;
      long long cnt = 1;
      ScalarMeanVar smv = ScalarMeanVar();
      ScalarMeanVar smv_var = ScalarMeanVar();
      for (int i = 0; i < N; i ++){
        ScalarMeanVar pre_smv = smv;
        smv.add(samples[i]);
        if (smv.getBiasedVar() > mid){
          smv_var.add(pre_smv.getBiasedVar());
          smv.reset();
          smv.add(samples[i]);
          cnt ++;
        }
      }

      if (cnt < cur_count){
        right = mid;
        smv_var.add(smv.getBiasedVar());
        my_vRAC = smv_var.getMax()/sigmas;
        my_vRWC = cnt*cnt*smv_var.getMax()/sigmas;
        my_cnt = cnt;
      } else left = mid;
    }
    double mid = (left+right)/2;
    exp.write("p", mid, my_cnt);
    exp.write("rac", mid, my_vRAC);
    exp.write("rwc", mid, my_vRWC);
    cur_count = my_cnt;
  }
}

VectorXd samples;
double sigmas;
double pmax, pmin;

double fv (double x[], int nvars){
  UNUSED(nvars);
  double beta = x[0];
  int cnt = 1;
  ScalarMeanVar smv = ScalarMeanVar();
  ScalarMeanVar smv_var = ScalarMeanVar();
  for (int i = 0; i < (int) samples.size(); i ++){
    ScalarMeanVar pre_smv = smv;
    smv.add(samples[i]);
    if (smv.getBiasedVar() > beta){
      smv_var.add(pre_smv.getBiasedVar());
      smv.reset();
      smv.add(samples[i]);
      cnt ++;
    }
  }
  if (cnt > pmax) return 24;
  if (cnt < pmin) return 24;
  smv_var.add(smv.getBiasedVar());
  double vRWC = cnt*cnt*smv_var.getMax()/sigmas;
  return vRWC;
}

double f (double beta){
  int cnt = 1;
  ScalarMeanVar smv = ScalarMeanVar();
  ScalarMeanVar smv_var = ScalarMeanVar();
  for (int i = 0; i < (int) samples.size(); i ++){
    ScalarMeanVar pre_smv = smv;
    smv.add(samples[i]);
    if (smv.getBiasedVar() > beta){
      smv_var.add(pre_smv.getBiasedVar());
      smv.reset();
      smv.add(samples[i]);
      cnt ++;
    }
  }
  if (cnt > pmax) return 24;
  if (cnt < pmin) return 24;
  smv_var.add(smv.getBiasedVar());
  cout << cnt << " " << smv_var.getMax() << endl;
  double vRWC = cnt*cnt*smv_var.getMax()/sigmas;
  return vRWC;
}

void dlv_iterative_test(){
  DetExp exp = DetExp("DLV3");
  int N = 100000;
  VectorXd samples;
  double mean = 0;
  double var = 100;
  double a = mean - sqrt(3*var);
  double b = mean + sqrt(3*var);
  UNUSED(a);
  UNUSED(b);
  Normal dist = Normal(mean, var);
  // Uniform dist = Uniform(a, b);
  dist.sample(samples, N);
  sort(samples.begin(), samples.end());
  ScalarMeanVar whole;
  for (int i = 0; i < (int) samples.size(); i ++){
    whole.add(samples[i]);
  }
  double sigmas = whole.getBiasedVar();
  double upper = 0.25*(whole.getMax() - whole.getMin())*(whole.getMax() - whole.getMin());
  double lower = 0;
  int cur_count = N;

  while (cur_count > 1){
    double left = lower;
    double right = upper;
    double my_cnt, my_vRAC, my_vRWC;
    while (abs(right-left)>1e-6){
      double mid = (left+right)/2;
      long long cnt = 1;
      ScalarMeanVar smv = ScalarMeanVar();
      ScalarMeanVar smv_var = ScalarMeanVar();
      for (int i = 0; i < N; i ++){
        ScalarMeanVar pre_smv = smv;
        smv.add(samples[i]);
        if (smv.getBiasedVar() > mid){
          smv_var.add(pre_smv.getBiasedVar());
          smv.reset();
          smv.add(samples[i]);
          cnt ++;
        }
      }

      if (cnt < cur_count){
        right = mid;
        smv_var.add(smv.getBiasedVar());
        my_vRAC = smv_var.getMax()/sigmas;
        my_vRWC = cnt*cnt*smv_var.getMean()/sigmas;
        my_cnt = cnt;
      } else left = mid;
    }

    priority_queue<pair<double, pair<int, int>>> pq;
    pq.push({0, {0, N-1}});

    double g_ratio = 0.1;

    while (pq.size() < my_cnt){
      auto p = pq.top().second;
      pq.pop();
      int l = p.first;
      int r = p.second;
      int n = r-l+1;
      ScalarMeanVar smv = ScalarMeanVar();
      for (int i = l; i <= r; i ++){
        smv.add(samples(i));
      }
      double beta = smv.getBiasedVar() * g_ratio * g_ratio * 13.5;
      smv.reset();
      int start = l;
      for (int i = l; i <= r; i ++){
        double total_var = smv.getM2();
        smv.add(samples(i));
        if (smv.getBiasedVar() > beta){
          pq.push({total_var, {start, i-1}});
          start = i;
          smv.reset();
          smv.add(samples(i));
        }
      }
      pq.push({smv.getM2(), {start, r}});
    }
    int this_cnt = pq.size();
    ScalarMeanVar smv_var = ScalarMeanVar();
    while (pq.size() > 0){
      auto p = pq.top().second;
      int n = p.second - p.first + 1;
      smv_var.add(pq.top().first/n);
      pq.pop();
    }

    my_vRAC = this_cnt*this_cnt*smv_var.getMean()/sigmas;
    double mid = (left+right)/2;
    exp.write("p", mid, my_cnt);
    exp.write("rac", mid, my_vRAC);
    exp.write("rwc", mid, my_vRWC);
    cur_count = my_cnt;
  }
}

void dlv_opt(){
  DetExp exp = DetExp("DLV_OPT");
  int N = 1000000;
  double group_ratio = 0.001;
  pmax = (int) ceil(N*group_ratio);
  pmin = (int) ceil(pmax / 13.5);
  double mean = 0;
  double var = 10000;
  double a = mean - sqrt(3*var);
  double b = mean + sqrt(3*var);
  UNUSED(a);
  UNUSED(b);
  Normal dist = Normal(mean, var);
  // Uniform dist = Uniform(a, b);
  dist.sample(samples, N);
  sort(samples.begin(), samples.end());
  ScalarMeanVar whole;
  for (int i = 0; i < (int) samples.size(); i ++){
    whole.add(samples[i]);
  }
  sigmas = whole.getBiasedVar();

  // double ibeta[1] = {1};
  // double ebeta[1];
  // double rho = 0.5;
  // double eps = 1e-8;
  // double itermax = 1000;
  // int iter = hooke(1, ibeta, ebeta, rho, eps, itermax, f);
  // cout << ebeta[0] << " " << iter << endl;

  double x1, x2;
  double dif = samples[N-1] - samples[0];
  double left = 0;
  double eps = 1e-5;
  double right = 0.25*dif*dif; // Popoviciu inequality
  // while (abs(left-right) > 1e-8){
  //   double mid1 = (2*left+right)/3;
  //   double mid2 = (left+2*right)/3;
  //   double fx1 = brent::local_min(0, mid1, eps, f, x1);
  //   double fx2 = brent::local_min(0, mid2, eps, f, x2);
  //   cout << left << " " << mid1 << " " << mid2 << " " << right << endl;
  //   cout << fx1 << " " << fx2 << endl;
  //   if (fx1 <= fx2) right = mid2;
  //   else left = mid1;
  // }
  // double mid = (left+right)/2;
  double x;
  double fx = brent::local_min(0, sigmas/(pmin*pmin), 1e-12, f, x);
  cout << x << " " << fx << endl;
}

int main(){
  // dlv_graph_test();
  dlv_iterative_test();
  // dlv_opt();
}