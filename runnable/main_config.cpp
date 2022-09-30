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

using namespace pb;

void main1(){
  // 19 20
  string file_name = "/home/alm818/model.lp";
  for (int seed = 1; seed <= 100; seed ++){
    vector<string> cols = {"tmass_prox", "j", "h", "k"};
    DetProb prob; prob.tableGenerate("ssds", cols, false, 200000, seed);
    prob.boundGenerate(1000, 0, 2);
    // vector<string> cols = {"price", "quantity", "discount", "tax"};
    // DetProb prob; prob.tableGenerate("tpch", cols, true, 10000, 5000, 0.5, seed);
    prob.normalizeObjective();
    cout << seed << "\n";
    // gs.writeModel(file_name);
    DualReducer dr = DualReducer(4, prob);
    printf("%.2f %.8f\n", dr.exe_lp, dr.lp_score);
    printf("%.2f %.8f\n", dr.exe_ilp, dr.ilp_score);
    GurobiSolver gs = GurobiSolver(prob);
    gs.solveIlp();
    DualReducer dr2 = DualReducer(4, prob, gs.ilp_sol);
    printf("%.2f %.8f\n", dr2.exe_ilp, dr2.ilp_score);
    if (dr.status != Found && dr2.status != Found && gs.ilp_status != Found) continue;
    double a, b, c;
    printf("%.2f %.8f %.4f%% %.4f%%\n", gs.getIlpTime(), gs.ilp_score, (dr.lp_score - dr.ilp_score)*100/(dr.lp_score - gs.ilp_score), (dr.lp_score - dr2.ilp_score)*100/(dr.lp_score - gs.ilp_score));
    a = fabs((dr.ilp_score - dr.lp_score) / dr.lp_score)*100;
    b = fabs((dr2.ilp_score - dr.lp_score) / dr.lp_score)*100;
    c = fabs((gs.ilp_score - dr.lp_score) / dr.lp_score)*100;
    printf("%.4f%% %.4f%% %.4f%% %.4f%% %.4f%%\n", a, b, c, a/c*100, b/c*100);
    Checker ch = Checker(prob, dr.kEpsilon);
    cout << solMessage(dr.status) << " " << solMessage(dr2.status) << endl;
    cout << feasMessage(ch.checkIlpFeasibility(dr.ilp_sol)) << " " << feasMessage(ch.checkIlpFeasibility(dr2.ilp_sol)) << endl;
    cout << ch.getScore(dr.ilp_sol) << " " << ch.getScore(dr2.ilp_sol) << " " << ch.getScore(gs.ilp_sol) << endl;
    if (dr.ilp_sol.size() < 100){
      cout << solCombination(dr.ilp_sol) << endl;
      cout << solCombination(dr2.ilp_sol) << endl;
      cout << solCombination(gs.ilp_sol) << endl;
    }
  }
}

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
        DynamicLowVariance dlv = DynamicLowVariance(group_ratio);
        dlv.partition(real_table_name, partition_name);
        dlv.pro.print();
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
  DynamicLowVariance dlv = DynamicLowVariance(group_ratio);
  dlv.partition(table_name, partition_name);
  dlv.pro.print();
}

void main6(){
  string table_name = "tpch_6_1";
  string partition_name = "P0";
  string obj_col = "price";
  bool is_maximize = true;
  vector<string> cols = {"quantity", "discount", "tax"};
  vector<int> consSense = {LowerBounded, UpperBounded, Bounded};
  LsrProb prob = LsrProb(table_name, partition_name, obj_col, is_maximize, cols, consSense);
  double E = 50;
  double alpha = 0;
  double hardness = 1;
  prob.boundGenerate(E, alpha, hardness);
  LayeredSketchRefine lsr = LayeredSketchRefine(kPCore, prob);
  lsr.pro.print();
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

int main(){
  // main4();
  //main5();
  // main5p5();
  main6();
  //main7();
}