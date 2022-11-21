#include "pb/det/det_sql.h"
#include "pb/det/det_exp.h"
#include "pb/det/det_prob.h"
#include "pb/det/lsr_prob.h"
#include "pb/det/dlv.h"
#include "pb/det/lsr.h"
#include "pb/det/sr.h"

#include "pb/core/checker.h"
#include "pb/core/gurobi_solver.h"
#include "pb/core/dual_reducer.h"

#include "pb/util/udeclare.h"
#include "pb/util/udebug.h"

using namespace pb;

int main() {
    KDTree kt;

    //vector<string> cols = {"tmass_prox", "j", "h", "k"};
    //kt.partitionTable("ssds_6_1", "s5000_d2", cols, 5000, 2);
    //vector<string> att_cols = {"j", "h", "k"};
    
    // vector<string> cols = {"price", "quantity", "discount", "tax"};
    // kt.partitionTable("tpch_6_1", "s5000_d100", cols, 5000, 100);
    vector<string> att_cols = {"quantity", "discount", "tax"};
    
    vector<int> att_senses = {0, 1, 2};
    DetSql det_sql = DetSql("tpch_6_1", "price", false, att_cols, att_senses, true, 1);
    //det_sql.addFilterWithRatio("j", 0.9, 0);
    //det_sql.addFilterWithRatio("k", 0.9, 2);
    double E = 150;
    double h = 5;
    LsrProb lsr_prob = LsrProb(det_sql, "s5000_d100", 1);
    lsr_prob.generateBounds(E, 0, h);
    fmt::print("E: {}\tHardness: {}\tis_maximizing: {}\n", E, h, lsr_prob.det_sql.is_maximize);
    DetProb det_prob = DetProb(det_sql, -1, 1);
    det_prob.copyBounds(lsr_prob.bl, lsr_prob.bu, lsr_prob.cl, lsr_prob.cu);
    
    SketchRefine sr = SketchRefine(lsr_prob);
    map<long long, long long> sol;
    sr.refine(sol);
    DualReducer dr = DualReducer(8, det_prob);

    //LsrChecker lsrCh = LsrChecker(lsr_prob);
    Checker ch = Checker(det_prob);
    double SR_score, DR_score, LP_score;
    SR_score = ch.getScore(sol);
    DR_score = ch.getScore(dr.ilp_sol);
    LP_score = ch.getScore(dr.lp_sol);

    fmt::print("DualReducer LP score: {}\n\n", LP_score);

    fmt::print("SketchRefine feasibility: {}\n", feasMessage(ch.checkIlpFeasibility(sol)));
    fmt::print("SketchRefine ILP score: {}\n", SR_score);
    fmt::print("SketchRefine pctError: {}\n", pctError(SR_score, LP_score));
    // GurobiSolver gs = GurobiSolver(det_prob, t       rue);
    // fmt::print("Query feasibility: {}\n",feasMessage() gs.hasIlpSolution());
    
    
    fmt::print("DualReducer feasibility: {}\n", feasMessage(ch.checkIlpFeasibility(dr.ilp_sol)));
    fmt::print("DualReducer ILP score: {}\n", DR_score);
    fmt::print("DualReducer pctError: {}\n", pctError(DR_score, LP_score));
}