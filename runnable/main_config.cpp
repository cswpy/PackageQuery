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
    vector<string> cols = {"tmass_prox", "j", "h", "k"};
    DetProb prob; prob.tableGenerate("ssds", cols, false, 200000, 1000, 0.8, seed);
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