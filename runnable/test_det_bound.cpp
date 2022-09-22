#include "pb/util/udeclare.h"
#include "pb/util/uconfig.h"
#include "pb/util/udebug.h"
#include "pb/det/det_prob.h"
#include "pb/core/dual_reducer.h"
#include "pb/core/gurobi_solver.h"
#include "pb/det/det_bound.h"
#include "pb/util/upostgres.h"

using namespace pb;

void test_hardness_distribution(){
  fstream out_file; out_file.open(kProjectHome + separator() + "plots" + separator() + "hardness.csv", std::ios::out);
  double E = 50;
  vector<string> table_names = {"ssds", "tpch"};
  vector<double> rhos = {0.1, 0.3, 0.5, 0.7, 0.9};
  vector<double> alphas = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};
  int repeats = 1000;
  PgManager pg = PgManager();
  for (string table_name : table_names){
    if (!pg.checkStats(table_name)){
      vector<string> cols = pg.getNumericCols(table_name);
      Stat *stat = pg.computeStats(table_name, cols);
      pg.writeStats(table_name, stat);
      delete stat;
    }
    Stat *stat = pg.readStats(table_name);
    vector<int> consSense = {LowerBounded, UpperBounded, Bounded};
    vector<double> means = {stat->mean(1), stat->mean(2), stat->mean(3)};
    vector<double> vars = {stat->getVar(1), stat->getVar(2), stat->getVar(3)};
    int m = 3;
    VectorXd bl (m); VectorXd bu (m);
    DetBound detBound = DetBound(consSense, means, vars);
    for (double rho : rhos){
      for (double alpha : alphas){
        for (int i = 0; i < repeats; i ++){
          double hard = detBound.sample(E, rho, alpha, bl, bu);
          string s = fmt::format("{},{:.2Lf},{:.2Lf},{:.8Lf}\n", table_name, rho, alpha, hard);
          cout << s;
          out_file << s;
          out_file.flush();
        }
      }
    }
    delete stat;
  }
  out_file.close();
}

void test_hardness(){
  fstream out_file; out_file.open(kProjectHome + separator() + "plots" + separator() + "dual_reducer_hardness.csv", std::ios::out);
  double E = 50;
  int N = 20000;
  vector<string> table_names = {"ssds", "tpch"};
  vector<double> alphas = {0.1, 0.5, 0.9};
  vector<double> hardnesses = {1, 3, 5, 7, 9, 11, 13};
  int repeats = 20;
  for (string table_name : table_names){
    DetProb prob;
    if (table_name == "ssds"){
      vector<string> cols = {"tmass_prox", "j", "h", "k"};
      prob.tableGenerate(table_name, cols, false, N);
    } else if (table_name == "tpch"){
      vector<string> cols = {"price", "quantity", "discount", "tax"};
      prob.tableGenerate(table_name, cols, true, N);
    }
    for (double hardness : hardnesses){
      for (double alpha : alphas){
        for (int i = 0; i < repeats; i ++){
          double rho = prob.boundGenerate(E, alpha, hardness);
          DualReducer dr = DualReducer(4, prob);
          GurobiSolver gs = GurobiSolver(prob);
          gs.solveLp();
          gs.solveIlp();
          string s = fmt::format("{},{:.2Lf},{:.2Lf},{},{},{},{:.8Lf},{:.8Lf},{:.8Lf},{:.2Lf},{:.2Lf},{:.2Lf}\n", 
            table_name, alpha, hardness,
            gs.lp_status == Found, dr.status == Found, gs.ilp_status == Found,
            gs.lp_score, dr.ilp_score, gs.ilp_score,
            gs.getLpTime(), dr.exe_ilp, gs.getIlpTime());
          cout << s;
          out_file << s;
          out_file.flush();
        }
      }
    }
  }
  out_file.close();
}

int main(){
  // test_hardness_distribution();
  test_hardness();
}