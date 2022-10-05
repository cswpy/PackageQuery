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

int main(){
  // main4();
  //main5();
  // main5p5();
  // main6();
  //main7();
  // main8();
  // main9();
  main10();
}