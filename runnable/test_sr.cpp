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

void test_ssds(string table_name, double size, double diameter, int order, DetExp &exp_size) {
    vector<int> att_senses = {0, 1, 2}; // Lower bounded, Upper bounded, Both bounded
    KDTree kt;
    vector<string> cols = {"tmass_prox", "j", "h", "k"};
    string partition_name = fmt::format("s{}_d{}", size, diameter);
    kt.partitionTable(table_name, partition_name, cols, size, (double)diameter);
    exp_size.write("ssds_kd_query", order, kt.exec_query);
    exp_size.write("ssds_kd_tree_building", order, kt.exec_tree_building);
    exp_size.write("ssds_kd_storage", order, kt.exec_storage);
    exp_size.write("ssds_kd_total", order, kt.exec_kd);

    vector<string> att_cols = {"j", "h", "k"};
    DetSql det_sql = DetSql(table_name, "tmass_prox", false, att_cols, att_senses, true, 1);

    DetExp exp = DetExp(table_name+"_"+partition_name);
    exp.o = order;
    exp.E = 50;
    vector<int> hardness = {1, 2, 3, 5, 7, 8, 10};
    exp.partition_name = partition_name;
    for(auto H : hardness) {
        exp.H = H;
        LsrProb lsr_prob = LsrProb(det_sql, partition_name, 1);
        lsr_prob.generateBounds(exp.E, 0, exp.H);
        SketchRefine sr = SketchRefine(lsr_prob);
        map<long long, long long> sol;
        bool status = sr.sketchAndRefine(sol);
        if (status) {
            exp.write(fmt::format("SR_{}", exp.E), exp.H, sr.exec_sr);
        }
    }

}

void test_ssds_5() {
    vector<int> att_senses = {0, 1, 2}; // Lower bounded, Upper bounded, Both bounded
    string partition_name = fmt::format("s{}_d{}", 500, 5);

    vector<string> att_cols = {"j", "h", "k"};
    DetSql det_sql = DetSql("ssds_5_1", "tmass_prox", false, att_cols, att_senses, true, 1);

    DetExp exp = DetExp("ssds_5_1_s500_d5_h");
    exp.o = 5;
    exp.E = 50;
    exp.partition_name = partition_name;
    int num_solvable_query = 0, num_solved_by_sr = 0;
    double total_time_elapsed = 0.0;
    for(auto H : exp.H8) {
        while(num_solvable_query<25) { // Varying seeds
            exp.H = H;
            exp.seed = num_solvable_query+1;
            LsrProb lsr_prob = LsrProb(det_sql, partition_name, exp.seed);
            lsr_prob.generateBounds(exp.E, 0, exp.H);
            DetProb det_prob = DetProb(det_sql, -1, exp.seed);
            det_prob.copyBounds(lsr_prob.bl, lsr_prob.bu, lsr_prob.cl, lsr_prob.cu);
            GurobiSolver gs = GurobiSolver(det_prob, true);
            gs.solveIlp();
            Checker ch = Checker(det_prob);
            if(ch.checkIlpFeasibility(gs.ilp_sol) != 1) continue;
            num_solvable_query++;
            SketchRefine sr = SketchRefine(lsr_prob);
            map<long long, long long> sol;
            bool status = sr.sketchAndRefine(sol);
            if (status) {
                num_solved_by_sr++;
                total_time_elapsed += sr.exec_sr;
            }
        }
        double avg_time_elapsed = total_time_elapsed / num_solved_by_sr;
        exp.write(fmt::format("SR_{}", exp.E), exp.H, avg_time_elapsed);
        exp.write("SR_accuracy", exp.H, num_solved_by_sr);
    }
}

void main2() {
    vector<int> att_senses = {0, 1, 2}; // Lower bounded, Upper bounded, Both bounded
    
    KDTree kt;

    // vector<string> cols = {"tmass_prox", "j", "h", "k"};
    // kt.partitionTable("ssds_5_1", "s500_d5", cols, 500, 5);
    // vector<string> att_cols = {"j", "h", "k"};
    // DetSql det_sql = DetSql("ssds_5_1", "tmass_prox", false, att_cols, att_senses, true, 1);

    vector<string> cols = {"price", "quantity", "discount", "tax"};
    kt.partitionTable("tpch_5_1", "s500", cols, 500, DBL_MAX);
    vector<string> att_cols = {"quantity", "discount", "tax"};
    DetSql det_sql = DetSql("tpch_5_1", "price", false, att_cols, att_senses, true, 1);
    
    //det_sql.addFilterWithRatio("j", 0.9, 0);
    //det_sql.addFilterWithRatio("k", 0.9, 2);

    double E = 50; // 50 or 500 for tpc-h
    double h = 7;
    LsrProb lsr_prob = LsrProb(det_sql, "s500", 1);
    lsr_prob.generateBounds(E, 0, h);
    fmt::print("E: {}\tHardness: {}\tis_maximizing: {}\n\n", E, h, lsr_prob.det_sql.is_maximize);
    DetProb det_prob = DetProb(det_sql, -1, 1);
    det_prob.copyBounds(lsr_prob.bl, lsr_prob.bu, lsr_prob.cl, lsr_prob.cu);
    
    SketchRefine sr = SketchRefine(lsr_prob);
    map<long long, long long> sol;
    bool status = sr.sketchAndRefine(sol);
    DualReducer dr = DualReducer(8, det_prob);

    //LsrChecker lsrCh = LsrChecker(lsr_prob);
    Checker ch = Checker(det_prob);

    double SR_score, DR_score, LP_score;
    
    DR_score = ch.getScore(dr.ilp_sol);
    LP_score = ch.getScore(dr.lp_sol);

    if(status) {
        SR_score = ch.getScore(sol);
        fmt::print("SketchRefine feasibility: {}\n", feasMessage(ch.checkIlpFeasibility(sol)));
        fmt::print("SketchRefine ILP score: {}\n", SR_score);
        fmt::print("SketchRefine pctError: {}\n\n", pctError(SR_score, LP_score));
    }

    GurobiSolver gs = GurobiSolver(det_prob, true);
    gs.solveIlp();
    fmt::print("Gurobi feasibility: {}\n", feasMessage(ch.checkIlpFeasibility(gs.ilp_sol)));
    fmt::print("Gurobi score: {}\n\n", gs.ilp_score);

    //fmt::print("DualReducer LP score: {}\n\n", LP_score);

    fmt::print("DualReducer feasibility: {}\n", feasMessage(ch.checkIlpFeasibility(dr.ilp_sol)));
    fmt::print("DualReducer ILP score: {}\n", DR_score);
    fmt::print("DualReducer pctError: {}\n\n", pctError(DR_score, LP_score));

    // GurobiSolver gs = GurobiSolver(det_prob, t       rue);
    // fmt::print("Query feasibility: {}\n",feasMessage() gs.hasIlpSolution());

}

int main() {
    // DetExp exp = DetExp("kd_partition");
    // vector<int> orders = {5, 6, 7, 8};
    // vector<int> size_reqs = {500, 5000, 50000, 500000};
    // vector<int> diameter_reqs = {3, 8, 15, 30};
    // for(int i=0; i<4; i++) {
    //     string table_name = fmt::format("ssds_{}_1", orders[i]);
    //     test_ssds(table_name, size_reqs[i], diameter_reqs[i], orders[i], exp);
    // }
    main2();
    //test_ssds_5();
}