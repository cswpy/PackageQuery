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

// void main1(){
//   // 19 20
//   string file_name = "/home/alm818/model.lp";
//   for (int seed = 1; seed <= 100; seed ++){
//     vector<string> cols = {"tmass_prox", "j", "h", "k"};
//     DetProb prob; prob.tableGenerate("ssds", cols, false, 200000, seed);
//     prob.boundGenerate(1000, 0, 2);
//     // vector<string> cols = {"price", "quantity", "discount", "tax"};
//     // DetProb prob; prob.tableGenerate("tpch", cols, true, 10000, 5000, 0.5, seed);
//     prob.normalizeObjective();
//     cout << seed << "\n";
//     // gs.writeModel(file_name);
//     DualReducer dr = DualReducer(4, prob);
//     printf("%.2f %.8f\n", dr.exe_lp, dr.lp_score);
//     printf("%.2f %.8f\n", dr.exe_ilp, dr.ilp_score);
//     GurobiSolver gs = GurobiSolver(prob);
//     gs.solveIlp();
//     DualReducer dr2 = DualReducer(4, prob, gs.ilp_sol);
//     printf("%.2f %.8f\n", dr2.exe_ilp, dr2.ilp_score);
//     if (dr.status != Found && dr2.status != Found && gs.ilp_status != Found) continue;
//     double a, b, c;
//     printf("%.2f %.8f %.4f%% %.4f%%\n", gs.getIlpTime(), gs.ilp_score, (dr.lp_score - dr.ilp_score)*100/(dr.lp_score - gs.ilp_score), (dr.lp_score - dr2.ilp_score)*100/(dr.lp_score - gs.ilp_score));
//     a = fabs((dr.ilp_score - dr.lp_score) / dr.lp_score)*100;
//     b = fabs((dr2.ilp_score - dr.lp_score) / dr.lp_score)*100;
//     c = fabs((gs.ilp_score - dr.lp_score) / dr.lp_score)*100;
//     printf("%.4f%% %.4f%% %.4f%% %.4f%% %.4f%%\n", a, b, c, a/c*100, b/c*100);
//     Checker ch = Checker(prob, dr.kEpsilon);
//     cout << solMessage(dr.status) << " " << solMessage(dr2.status) << endl;
//     cout << feasMessage(ch.checkIlpFeasibility(dr.ilp_sol)) << " " << feasMessage(ch.checkIlpFeasibility(dr2.ilp_sol)) << endl;
//     cout << ch.getScore(dr.ilp_sol) << " " << ch.getScore(dr2.ilp_sol) << " " << ch.getScore(gs.ilp_sol) << endl;
//     if (dr.ilp_sol.size() < 100){
//       cout << solCombination(dr.ilp_sol) << endl;
//       cout << solCombination(dr2.ilp_sol) << endl;
//       cout << solCombination(gs.ilp_sol) << endl;
//     }
//   }
// }

void main2(){
  vector<int> consSense = {UpperBounded, LowerBounded, Bounded, LowerBounded, UpperBounded};
  vector<double> means = {10, 20, 30, 40, 50};
  vector<double> vars = {100, 100, 150, 400, 100};
  DetBound detBound = DetBound(consSense, means, vars);
  double E = 50;
  double rho = 0.8;
  double alpha = 0.5;
  int m = consSense.size();
  VectorXd bl (m); bl.fill(-DBL_MAX);
  VectorXd bu (m); bu.fill(DBL_MAX);
  double hard = detBound.sample(E, rho, alpha, bl, bu);
  double minE;
  double minHard = detBound.minHardness(minE, bl, bu);
  cout << hard << endl;
  cout << minE << " " << minHard << endl;
  double this_rho = detBound.sampleHardness(E, alpha, 10, bl, bu);
  double my_hard = detBound.measureHardness(E, bl, bu);
  cout << this_rho << " " << my_hard << endl;
}

void main3(){
  string table_name = "ssds";
  PgManager pg = PgManager();
  pg.readStats(table_name);
}

void main4(){
  Synthetic syn = Synthetic();
  // long long N = 1000000000;
  // vector<string> cols = {"tmass_prox", "j", "h", "k"};
  // syn.createSubtable("ssds", 8, cols, 1);
  vector<string> cols = {"price", "quantity", "discount", "tax"};
  syn.createSubtable("tpch", 8, cols, 1);
  syn.pro.print();
}

void main5(){
  double group_ratio = 0.01;
  string partition_name = "P0";
  vector<string> table_names = {"tpch", "ssds"};
  vector<vector<int>> orders = {{6,7,8,9}, {6,7,8}};
  int seeds = 1;
  for (int t = 0; t < (int) table_names.size(); t ++){
    string table_name = table_names[t];
    for (int order : orders[t]){
      for (int seed = 1; seed <= seeds; seed++){
        string real_table_name = fmt::format("{}_{}_{}", table_name, order, seed);
        cout << table_name << " on order " << order << endl;
        DynamicLowVariance dlv = DynamicLowVariance(kPCore, group_ratio);
        dlv.partition(real_table_name, partition_name);
      }
    }
  }
  // dlv.dropTempTables();
  // dlv.dropPartition(table_name, partition_name);
}

void main5p5(){
  double group_ratio = 0.01;
  string partition_name = "P1";
  string table_name = "ssds_7_1";
  DynamicLowVariance dlv = DynamicLowVariance(kPCore, group_ratio);
  dlv.partition(table_name, partition_name);
}

void main6(){
  Profiler pro = Profiler({"DetProb"});
  for (int rep = 0; rep < 1; rep ++){
    int order = 9;
    string table_name = fmt::format("tpch_{}_1", order);
    string partition_name = "P0";
    string obj_col = "price";
    bool is_maximize = true;
    vector<string> att_cols = {"quantity", "discount", "tax"};
    vector<int> att_senses = {LowerBounded, UpperBounded, Bounded};
    int att_count = (int) att_cols.size();
    DetSql det_sql = DetSql(table_name, obj_col, is_maximize, att_cols, att_senses);
    LsrProb prob = LsrProb(det_sql, partition_name);
    double E = 50;
    double alpha = 0;
    double hardness = 8;
    prob.generateBounds(E, alpha, hardness);
    LayeredSketchRefine lsr = LayeredSketchRefine(kPCore, prob);
    #if DEBUG
      lsr.pro.print();
    #endif
    if (order <= 7){
      long long n = (long long) (2 * pow(10.0, order));
      pro.clock(0);
      DetProb det_prob = DetProb(det_sql, n);
      pro.stop(0);
      pro.print();
      memcpy(&det_prob.bl(0), &prob.bl(0), prob.bl.size()*sizeof(double));
      memcpy(&det_prob.bu(0), &prob.bu(0), prob.bu.size()*sizeof(double));
      det_prob.bl(att_count) = prob.cl;
      det_prob.bu(att_count) = prob.cu;
      det_prob.truncate();
      Checker ch = Checker(det_prob);
      DualReducer dr = DualReducer(kPCore, det_prob);
      double ground = dr.lp_score;
      cout << feasMessage(ch.checkLpFeasibility(lsr.lp_sol)) << " " << feasMessage(ch.checkIlpFeasibility(lsr.ilp_sol)) << " " << feasMessage(ch.checkLpFeasibility(dr.lp_sol)) << " " << feasMessage(ch.checkIlpFeasibility(lsr.ilp_sol)) << endl;
      cout << "LSR" << endl;
      cout << fmt::format("{:.8Lf} {:.8Lf} {:.4Lf}% {:.4Lf}%\n", lsr.lp_score, lsr.ilp_score, pctError(lsr.lp_score, ground), pctError(lsr.ilp_score, ground));
      cout << fmt::format("{:.2Lf}ms {:.2Lf}ms\n", lsr.exe_lp, lsr.exe_ilp);
      cout << "DR" << endl;
      cout << fmt::format("{:.8Lf} {:.8Lf} {:.4Lf}%\n", dr.lp_score, dr.ilp_score, pctError(dr.ilp_score, ground));
      cout << fmt::format("{:.2Lf}ms {:.2Lf}ms\n", dr.exe_lp, dr.exe_ilp);
    } else{
      cout << "LSR" << endl;
      LsrChecker ch = LsrChecker(prob);
      double ground = lsr.lp_score;
      cout << feasMessage(ch.checkLpFeasibility(lsr.lp_sol)) << " " << feasMessage(ch.checkIlpFeasibility(lsr.ilp_sol)) << endl;
      cout << fmt::format("{:.8Lf} {:.8Lf} {:.4Lf}% {:.4Lf}%\n", lsr.lp_score, lsr.ilp_score, pctError(lsr.lp_score, ground), pctError(lsr.ilp_score, ground));
      cout << fmt::format("{:.2Lf}ms {:.2Lf}ms\n", lsr.exe_lp, lsr.exe_ilp);
    }
    // print(lsr.lp_sol);
    // print(lsr.ilp_sol);
  }
}

void main7(){
  int m = 3;
  int n = 2;
  DetProb prob = DetProb(m, n);
  for (int i = 0; i < m; i ++){
    for (int j = 0; j < n; j ++){
      prob.A(i, j) = i*n+j;
    }
  }
  prob.bl(0) = 1;
  prob.bu(1) = 2;
  prob.truncate();
  cout << prob.A << endl;
  print(prob.bu);
  print(prob.bl);
  DetProb prob2 = prob;
  prob.A(0, 0) = 100;
  cout << prob2.A << endl;
}

void main8(){
  int n = 4;
  VectorXi a (n);
  for (int i = 0; i < n; i ++) a(i) = i+1;
  VectorXi b (n*2); b.fill(0);
  memcpy(&b(2), &a(0), n*sizeof(int));
  print(a);
  print(b);
}

void main9(){
  int r = 3;
  int c = 4;
  RMatrixXd A (r, c);
  MatrixXd B (r, c);
  for (int i = 0; i < r; i ++){
    for (int j = 0; j < c; j ++){
      A(i, j) = i * c + j + 1;
      B(i, j) = A(i, j);
    }
  }
  cout << A << endl << endl;
  VectorXd C(5);
  // Conclusion
  // MatrixXd is column major meaning increasing address corresponding to row-first then to column
  // RMatrixXd is row major meaning increasing address corresponding to column-first then to row
  memcpy(&C(0), &A(1,0), 3*sizeof(double));
  // If source is RMatrixXd m x sub_n
  // Dest is RMatrixXd m x n 
  // WE ARE GOOD HERE.
  print(C);
  // cout << B << endl;
}

bool is_finished = false;
double max_RAM = 0;
const int kSleepPeriod = 10; // In ms

void inspectRAM(){
  while (!is_finished){
    max_RAM = max(max_RAM, currentRAM());
    std::this_thread::sleep_for(std::chrono::milliseconds(kSleepPeriod));
  }
}

void main10(){
  // std::thread ram (inspectRAM);
  double group_ratio = 0.01;
  string partition_name = "P1";
  string table_name = "tpch_9_1";
  DynamicLowVariance dlv = DynamicLowVariance(kPCore, group_ratio);
  // dlv.dropAllPartitions();
  // dlv.dropTempTables();
  // dlv.partition(table_name, partition_name);
  // dlv.pro.print();
  // cout << dlv.exe << endl;
  // is_finished = true;
  // ram.join();
  // cout << "RAM:" << max_RAM << endl;
}

// void main11(){
//   double mu = 10;
//   double variance = 100;
//   Profiler pro = Profiler({"1", "2"});
//   Dist* normal = new Normal(mu, variance);
//   int N = 1000000;
//   pro.clock(0);
//   for (int i = 0; i < N-1; i ++){
//     normal->quantile(1.0/N*(i+1));
//     normal->cdf(i-N/2);
//   }
//   pro.stop(0);
//   pro.clock(1);
//   for (int i = 0; i < N-1; i ++){
//     normalQuantile(mu, variance, 1.0/N*(i+1));
//     normalCdf(mu, variance, i-N/2);
//   }
//   pro.stop(1);
//   // VectorXd x;
//   // normal->sample(x, 10, -1);
//   // print(x);
//   pro.print();
// }

// void main12(){
//   double mu = 0;
//   double sigma = 1;
//   Dist* dist = new LogNormal(mu, sigma);
//   VectorXd x;
//   dist->sample(x, 10, -1);
//   print(x);
//   cout << dist->left_cvar(1.0) << " " << dist->right_cvar(1.0) << endl;
//   cout << dist->left_cvar(0.1) << " " << dist->right_cvar(0.1) << endl; 
//   cout << dist->mean() << " " << dist->variance() << endl;
// }

void main13(){
  DetExp exp = DetExp("L2");
  exp.o = 7;
  exp.partition_name = "P3";
  DetSql det_sql = exp.generate();
  LsrProb lsr_prob = LsrProb(det_sql, exp.partition_name, exp.seed);
  lsr_prob.generateBounds(exp.E, exp.a, exp.H);
  LayeredSketchRefine lsr = LayeredSketchRefine(exp.C, lsr_prob);
  cout << solMessage(lsr.status) << " " << lsr.exe_ilp << endl;
}

vector<int> getBuckets(int from, int to, int bucket){
  assert(to > from);
  vector<int> res;
  int dif = to - from;
  if (dif <= bucket){
    res.resize(dif + 1);
    iota(res.begin(), res.end(), from);
    return res;
  }
  res.resize(bucket + 1);
  res[0] = from;
  div_t q = div(dif, bucket);
  for (int i = 0; i < q.rem; i ++) res[i+1] = res[i] + q.quot + 1;
  for (int i = q.rem; i <= bucket; i ++) res[i+1] = res[i] + q.quot;
  return res;
}

double my_quantile(map<int, double> &qtree, int N, double p){
  double bucket_width = 1/(N+1.0);
  double left_most = bucket_width;
  if (p <= left_most) return qtree[0]-(qtree[1] - qtree[0])*(left_most-p)/bucket_width;
  double right_most = N*bucket_width;
  if (p >= right_most) return qtree[N-1]+(qtree[N-1] - qtree[N-2])*(p-right_most)/bucket_width;
  int low = (int) floor(p / bucket_width);
  auto ptr = qtree.upper_bound(low);
  // cout << p << " " << bucket_width << " " << low << " " << ptr->first << endl;
  double sv = ptr->second;
  double s = (ptr->first+1)*bucket_width;
  double f = (prev(ptr)->first+1)*bucket_width;
  double fv = prev(ptr)->second;
  return fv + (sv-fv)*(p-f)/(s-f);
}

void main14(){
  double eps = 0.01;
  RandomQuantile rq (eps);
  LogNormal dist = LogNormal(0, 3);
  VectorXd samples;
  int N = 1000000;
  dist.sample(samples, N);
  for (auto v : samples) rq.feed(v);
  rq.finalize();
  double left = 0.01;
  double right = 0.01;
  double deps = 0.000001;
  double left_error = 0;
  double right_error = 0;
  int cnt = 0;
  double veps = 1e-4;
  for (double p = veps; p <= left; p += deps){
    left_error += abs(dist.quantile(p)-rq.query_for_value(p));
    cnt ++;
  }
  left_error /= cnt;
  cnt = 0;
  for(double p = 1-right; p <= 1-veps; p += deps){
    right_error += abs(dist.quantile(p)-rq.query_for_value(p));
    // cout << right_error << endl;
    cnt ++;
  }
  right_error /= cnt;
  cout << left_error << " " << right_error << endl;

  int estimated = (int) ceil(1/eps*pow(log2(1/eps), 1.5));
  sort(samples.begin(), samples.end());
  int bucket = 10;
  int limit = estimated;
  cout << estimated << endl;
  vector<int> buckets = getBuckets(0, N-1, bucket);
  vector<tuple<double, int, int>> pq (limit + bucket);
  for (int i = 0; i < bucket; i ++){
    double v = samples[buckets[i+1]] - samples[buckets[i]];
    pq[i] = {v, buckets[i], buckets[i+1]};
  }
  int heap_size = bucket;
  make_heap(pq.begin(), pq.begin()+heap_size);
  while (heap_size < limit){
    auto [v, from, to] = pq[0];
    pop_heap(pq.begin(), pq.begin()+heap_size);
    heap_size --;
    if (to - from <= 1) continue;
    buckets = getBuckets(from, to, bucket);
    for (int i = 0; i < bucket; i ++){
      double v = samples[buckets[i+1]] - samples[buckets[i]];
      pq[heap_size] = {v, buckets[i], buckets[i+1]};
      heap_size ++;
      push_heap(pq.begin(), pq.begin()+heap_size);
    }
  }
  map<int, double> qtree;
  for (int i = 0; i < heap_size; i ++){
    auto [v, from, to] = pq[i];
    qtree[from] = samples[from];
    qtree[to] = samples[to];
  }
  qtree[1] = samples[1];
  qtree[N-2] = samples[N-2];

  left_error = 0;
  right_error = 0;
  
  cnt = 0;
  for (double p = veps; p <= left; p += deps){
    left_error += abs(dist.quantile(p)-my_quantile(qtree, N, p));
    cnt ++;
  }
  left_error /= cnt;
  cnt = 0;
  for(double p = 1-right; p <= 1-veps; p += deps){
    // cout << p << " " << dist.quantile(p) << " " << my_quantile(qtree, N, p) << " " << rq.query_for_value(p) << endl;
    right_error += abs(dist.quantile(p)-my_quantile(qtree, N, p));
    // cout << right_error << endl;
    cnt ++;
  }
  right_error /= cnt;
  cout << left_error << " " << right_error << endl;  
  cout << cnt << " " << qtree[0] << " " << qtree[1] << " " << qtree[N-2] << " " << qtree[N-1] << endl;
}

void main15(){
  int n = 50;
  int seed = -1;
  if (seed < 0){
    random_device rd;
    seed = rd();
  }
  vector<VG*> dists (n);
  vector<double> ms (n);
  default_random_engine gen (seed);
  double my_mean = 0;
  double my_variance = 0;
  for (int i = 0; i < n; i ++){
    ms[i] = r8_uniform_ab(0, 1, seed);
    int choice = i4_uniform_ab(0, 0, seed);
    if (choice == 0){
      uniform_real_distribution dist1 (0.0, 10.0);
      uniform_real_distribution dist2 (50.0, 100.0);
      dists[i] = new Uniform(dist1(gen), dist2(gen));
    } else if (choice == 1){
      uniform_real_distribution dist1 (-10.0, 10.0);
      uniform_real_distribution dist2 (100.0, 10000.0);
      dists[i] = new Normal(dist1(gen), dist2(gen));
    } else{
      uniform_real_distribution dist1 (0.0, 0.25);
      uniform_real_distribution dist2 (0.0, 2.0);
      dists[i] = new LogNormal(dist1(gen), dist2(gen));
    }
  }
  vector<double> norm_ms = linearCombination(ms);
  for (int i = 0; i < n; i ++){
    Dist* o = (Dist*) dists[i];
    my_mean += norm_ms[i]*o->mean();
    my_variance += norm_ms[i]*norm_ms[i]*o->variance();
  }

  LinearVG lvg = LinearVG(dists, ms);
  MixtureVG mvg = MixtureVG(dists, ms);

  int R = 1000000;
  // print(lvg.ms);
  // print(mvg.ms);
  cout << my_mean << " " << my_variance << endl;
  cout << lvg.mean(R) << " " << lvg.variance(R) << endl;
  cout << mvg.mean(R) << " " << mvg.variance(R) << endl;
  // double d1 = wassersteinDistance(lvg, mvg, R, seed);
  // cout << d1 << endl;
  // Normal en = Normal(my_mean, my_variance);
  // double d2 = wassersteinDistance(lvg, en, R, seed);
  // cout << d2 << endl;
  for (int i = 0; i < n; i ++){
    delete dists[i];
  }
}

int main(){
  // main4();
  //main5();
  // main5p5();
  // main6();
  //main7();
  // main8();
  // main9();
  // main10();
  // main11();
  // main12();
  // main13();
  main14();
  // main15();
}