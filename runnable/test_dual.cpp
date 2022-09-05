#include "pb/util/udeclare.h"
#include "pb/util/udebug.h"
#include "pb/util/uconfig.h"
#include "pb/core/dual.h"
#include "pb/core/gurobi_solver.h"

using namespace pb;

void testAcc(){
  int N = 100000;
  int E = 50;
  double var = 10;
  double outlier_prob = 2.0;
  int T = 1000;
  for (int i = 0; i < T; i ++){
    DetProb prob; prob.uniformGenerate(N, E, var, outlier_prob);
    Dual dual = Dual(kPCore, prob);
    GurobiSolver gs = GurobiSolver(prob);
    gs.solveLp();
    fmt::print("{:.10Lf} {:.10Lf}\n", dual.score, gs.lp_score);
    assert(gs.lp_score - dual.score < 1e-6);
  }
}

void testRunningTime(){
  vector<int> Es = {10, 10, 1000, 1000};
  vector<double> vars = {1, 10000, 1, 10000};
  vector<int> Ns = {1000, 10000, 100000, 1000000, 10000000, 100000000};
  // vector<int> Ns = {100000000};
  vector<int> seeds = {1, 2, 3, 4};
  double outlier_prob = 2.0;
  for (int N : Ns){
    for (int q = 0 ; q < 8; q ++){
      int index = q % 4;
      for (int seed : seeds){
        DetProb prob;
        if (q < 4) prob.uniformGenerate(N, Es[index], vars[index], outlier_prob, true, false, false, seed);
        else prob.normalGenerate(N, Es[index], vars[index], outlier_prob, true, seed);
        Dual d1 = Dual(1, prob);
        Dual d4 = Dual(4, prob);
        Dual d80 = Dual(80, prob);
        string s;
        if (N <= 10000000){
          GurobiSolver gs = GurobiSolver(prob);
          gs.solveLp();
          s = fmt::format("{},{},{},{},{:.2Lf},{:.2Lf},{:.2Lf},{:.2Lf},{:.8Lf},{:.8Lf},{:.8Lf},{:.8Lf}\n", N, Es[index], vars[index], seed, d1.exe_solve, d4.exe_solve, d80.exe_solve, gs.getLpTime(), d1.score, d4.score, d80.score, gs.lp_score);
        } else{
          s = fmt::format("{},{},{},{},{:.2Lf},{:.2Lf},{:.2Lf},{:.2Lf},{:.8Lf},{:.8Lf},{:.8Lf},{:.8Lf}\n", N, Es[index], vars[index], seed, d1.exe_solve, d4.exe_solve, d80.exe_solve, -1.0, d1.score, d4.score, d80.score, -1.0);
        }
        cout << s;
      }
    }
  }
}

int main(){
  //testAcc();
  testRunningTime();
}