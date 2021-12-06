#include "pb/util/udeclare.h"
#include "pb/util/udebug.h"
#include "pb/util/uconfig.h"
#include "pb/core/dual.h"
#include "pb/core/gurobi_solver.h"

using namespace pb;

void testAcc(){
  int N = 10000;
  int E = 50;
  double var = 10000;
  int T = 1000;
  for (int i = 0; i < T; i ++){
    DetProb prob; prob.uniformGenerate(N, E, var);
    Dual dual = Dual(kPCore, prob);
    GurobiSolver gs = GurobiSolver(prob);
    gs.solveLp();
    fmt::print("{:.10Lf} {:.10Lf}\n", dual.score, gs.lp_score);
    assert(gs.lp_score - dual.score < 1e-6);
  }
}

int main(){
}