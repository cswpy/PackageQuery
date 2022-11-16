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
#include "pb/core/vg.h"

#include "pb/lib/random_quantile.h"
#include "pb/lib/common.h"

using namespace pb;

struct hash_pair {
    template <class T1, class T2>
    size_t operator()(const pair<T1, T2>& p) const
    {
        auto hash1 = hash<T1>{}(p.first);
        auto hash2 = hash<T2>{}(p.second);
 
        if (hash1 != hash2) {
            return hash1 ^ hash2;             
        }
         
        // If hash1 == hash2, their XOR is zero.
          return hash1;
    }
};

void main1(){
  int N = 1000000;
  double group_ratio = 0.01;
  int g = (int) ceil(N*group_ratio);
  VectorXd samples;
  double mean = 0;
  double var = 10000;
  double a = mean - sqrt(3*var);
  double b = mean + sqrt(3*var);
  Normal dist = Normal(mean, var);
  // Uniform dist = Uniform(a, b);
  dist.sample(samples, N);
  sort(samples.begin(), samples.end());
  set<tuple<double, double, double>> fuse_set;
  map<double, pair<ScalarMeanVar, int>> partition_map;
  unordered_map<pair<double, double>, double, hash_pair> fuse_key_map;
  for (int i = 0; i < N; i ++){
    ScalarMeanVar smv = ScalarMeanVar();
    smv.add(samples[i]);
    partition_map[i] = {smv, i};
    if (i < N-1){
      smv.add(samples[i+1]);
      double v = smv.getBiasedVar();
      fuse_set.emplace(v, i, i+1);
      fuse_key_map[{i, i+1}] = v;
    }
  }
  for (int i = 0; i < N-g; i ++){
    auto &[_, i1, i2] = *fuse_set.begin();
    // cout << "IT " << i << " " << _ << " " << i1 << " " << i2 << endl;
    auto it1 = partition_map.find(i1);
    auto &[smv1, si1] = it1->second;
    auto it2 = next(it1);
    auto &[smv2, si2] = it2->second;
    ScalarMeanVar smv = smv1;
    smv.add(smv2);
    double mid = (i1+i2)/2;
    // cout << "WTF 1 "<< fuse_set.size() << " " << partition_map.size() << endl;
    if (it1 != partition_map.begin()){
      auto pit1 = prev(it1);
      double pi1 = pit1->first;
      auto psmv1 = pit1->second.first;
      // if (psmv1.sample_count > N) cout << "HERE " << psmv1.sample_count << endl;
      ////
      // double oracle = -1;
      // for (auto &[_1,_2,_3] : fuse_set){
      //   if (_2 == pi1 && _3 == i1){
      //     oracle = _1;
      //     break;
      //   }
      // }
      // cout << oracle << " " << psmv1.getBiasedVar() << " " << pi1 << " " << i1 << endl;
      // double a = roundLf(oracle, precision);
      // double b = roundLf(psmv1.getBiasedVar(), precision);
      // cout << fmt::format("{:.20Lf} {:.20Lf} {}\n", a, b, a==b);
      ////
      fuse_set.erase({fuse_key_map[{pi1, i1}], pi1, i1});
      fuse_key_map.erase({pi1, i1});
      psmv1.add(smv);
      double v = psmv1.getBiasedVar();
      // if (isnan(v)) cout << "WTF " << psmv1.sample_count << endl;
      fuse_set.emplace(v, pi1, mid);
      fuse_key_map[{pi1, mid}] = v;
    }
    // cout << "WTF 2 "<< fuse_set.size() << " " << partition_map.size() << endl;

    auto pit2 = next(it2);
    if (pit2 != partition_map.end()){
      double pi2 = pit2->first;
      auto psmv2 = pit2->second.first;
      // if (psmv2.sample_count > N) cout << "HERE " << psmv2.sample_count << endl;
      ////
      // double oracle = -1;
      // for (auto &[_1,_2,_3] : fuse_set){
      //   if (_2 == i2 && _3 == pi2){
      //     oracle = _1;
      //     break;
      //   }
      // }
      // cout << oracle << " " << psmv2.getBiasedVar() << " " << i2 << " " << pi2 << endl;
      // double a = roundLf(oracle, precision);
      // double b = roundLf(psmv2.getBiasedVar(), precision);
      // cout << fmt::format("{:.20Lf} {:.20Lf} {}\n", a, b, a==b);
      ////
      fuse_set.erase({fuse_key_map[{i2, pi2}], i2, pi2});
      fuse_key_map.erase({i2, pi2});
      psmv2.add(smv);
      double v = psmv2.getBiasedVar();
      // if (isnan(v)) cout << "WTF " << psmv2.sample_count << endl;
      fuse_set.emplace(v, mid, pi2);
      fuse_key_map[{mid, pi2}] = v;
    }
    // cout << "WTF 3 "<< fuse_set.size() << " " << partition_map.size() << endl;
    fuse_set.erase({_, i1, i2});
    fuse_key_map.erase({i1, i2});
    partition_map[mid] = {smv, si1};
    partition_map.erase(i1);
    partition_map.erase(i2);
    // cout << "WTF 4 "<< fuse_set.size() << " " << partition_map.size() << endl;
    // if (fuse_set.size() + 1 != partition_map.size()){
    //   break;
    // }
  }

  ScalarMeanVar smv_var = ScalarMeanVar();

  vector<double> delimiters;
  vector<int> sizes;
  delimiters.reserve(partition_map.size()-1);
  sizes.reserve(partition_map.size());
  sizes.push_back(partition_map.begin()->second.first.sample_count);
  smv_var.add(partition_map.begin()->second.first.getBiasedVar());
  for (auto it = next(partition_map.begin()); it != partition_map.end(); it ++){
    int index = it->second.second;
    double d = (samples[index] + samples[index-1]) / 2.0;
    delimiters.push_back(d);
    sizes.push_back(it->second.first.sample_count);
    smv_var.add(it->second.first.getBiasedVar());
  }
  // print(sizes);

  double left = 0;
  double dif = samples[N-1] - samples[0];
  double right = 0.25*dif*dif; // Popoviciu inequality
  while (abs(left-right) > 1e-8){
    double mid = (left+right)/2;
    int cnt = 1;
    ScalarMeanVar smv = ScalarMeanVar();
    for (int i = 0; i < N; i ++){
      smv.add(samples[i]);
      if (smv.getBiasedVar() > mid){
        smv.reset();
        smv.add(samples[i]);
        cnt ++;
      }
    }
    // cout << cnt << " " << g << " " << left << " " << mid << " " << right << endl;
    if (cnt < g) right = mid;
    else left = mid;
  }
  ScalarMeanVar smv_opt_var = ScalarMeanVar();
  double opt_var = (left+right) / 2;
  vector<int> opt_sizes (g, 0);
  int index = 0;
  ScalarMeanVar smv = ScalarMeanVar();
  for (int i = 0; i < N; i ++){
    ScalarMeanVar tmp_smv = smv;
    smv.add(samples[i]);
    if (smv.getBiasedVar() > opt_var){
      smv_opt_var.add(tmp_smv.getBiasedVar());
      smv.reset();
      smv.add(samples[i]);
      index ++;
    }
    opt_sizes[index] ++;
  }
  smv_opt_var.add(smv.getBiasedVar());
  cout << smv_var.getMean() << " " << smv_var.getMax() << endl;
  cout << smv_opt_var.getMean() << " " << smv_opt_var.getMax() << endl;
  // print(opt_sizes);
}

int main(){
  main1();
}