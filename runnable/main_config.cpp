#include "pb/util/udeclare.h"
#include "pb/util/uconfig.h"
#include "pb/util/udebug.h"
#include "pb/det/det_prob.h"
#include "pb/core/dual_reducer.h"
#include "pb/core/gurobi_solver.h"

using namespace pb;

int main(){
  // 19 20
  string file_name = "/home/alm818/model.lp";
  for (int seed = 1; seed <= 100; seed ++){
    // vector<string> cols = {"tmass_prox", "j", "h", "k"};
    // DetProb prob; prob.tableGenerate("ssds", cols, false, 100000, 10, 0.5, seed);
    vector<string> cols = {"price", "quantity", "discount", "tax"};
    DetProb prob; prob.tableGenerate("tpch", cols, true, 10000, 10, 0.01, seed);
    // prob.normalizeObjective();
    cout << seed << "\n";
    GurobiSolver gs = GurobiSolver(prob);
    gs.solveIlp();
    // gs.writeModel(file_name);
    DualReducer dr = DualReducer(4, prob);
    DualReducer dr2 = DualReducer(4, prob, gs.ilp_sol);
    if (dr.status != Found && dr2.status != Found && gs.ilp_status != Found) continue;
    printf("%.2f %.8f\n", dr.exe_lp, dr.lp_score);
    printf("%.2f %.8f\n", dr.exe_ilp, dr.ilp_score);
    printf("%.2f %.8f\n", dr2.exe_ilp, dr2.ilp_score);
    printf("%.2f %.8f %.4f%% %.4f%%\n", gs.getIlpTime(), gs.ilp_score, (dr.lp_score - dr.ilp_score)*100/(dr.lp_score - gs.ilp_score), (dr.lp_score - dr2.ilp_score)*100/(dr.lp_score - gs.ilp_score));
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