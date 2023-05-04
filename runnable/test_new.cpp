#include "pb/util/udeclare.h"
#include "pb/det/kd_tree.h"
#include "pb/det/det_sql.h"
#include "pb/det/det_prob.h"
#include "pb/core/gurobi_solver.h"
#include "pb/util/udebug.h"
#include "pb/det/incremental_eval.h"

using namespace pb;

void print_sol(VectorXd a) {
    for(int i=0; i<a.size(); i++) {
        if (a(i)!=0) {
            fmt::print("({}), {}  ", i, a(i));
        }
    }
    fmt::print("\n");
}

void print_sol(map<long long, long long> sol) {
    for(auto it=sol.begin(); it!=sol.end(); it++){
        fmt::print("({}), {}  ", it->first, it->second);
    }
    fmt::print("\n");
}

int test_single_ILP() {
    vector<int> att_senses = {0, 1, 2}; // Lower bounded, Upper bounded, Both bounded

    // vector<string> cols = {"tmass_prox", "j", "h", "k"};
    // //kt.partitionTable("ssds_3_1", "s500_d5", cols, 500, 5);
    // vector<string> att_cols = {"j", "h", "k"};
    // DetSql det_sql = DetSql("ssds_5_1", "tmass_prox", false, att_cols, att_senses, true, 1);

    vector<string> cols = {"price", "quantity", "discount", "tax"};
    vector<string> att_cols = {"quantity", "discount", "tax"};
    DetSql det_sql = DetSql("tpch_5_1", "price", false, att_cols, att_senses, true, 1);
    
    //det_sql.addFilterWithRatio("j", 0.9, 0);
    //det_sql.addFilterWithRatio("k", 0.9, 2);

    // IOFormat printVecFmt(2, 0, ",", " ", "", "", "[", "]");

    double E = 10; // 50 or 500 for tpc-h
    double h = 1;
    // LsrProb lsr_prob = LsrProb(det_sql, "s500", 1);
    DetProb q1 = DetProb(det_sql, -1, 1);
    q1.generateBounds(E, 0, h);
    DetProb q2 = DetProb(det_sql, -1, 1);
    q2.generateBounds(E, 0, h+0.1);
    q1.display();
    q2.display();

    GurobiSolver gs = GurobiSolver(q1);
    gs.solveLp();
    gs.solveIlp(0, -1);
    GurobiSolver gs2 = GurobiSolver(q2);
    gs2.solveLp();
    gs2.solveIlp(0, -1);

    fmt::print("Q1 Gurobi LP status: {}\t Score: {}\n", feasMessage(gs.checkLpFeasibility(gs.lp_sol)), gs.lp_score);
    for (int i=0; i<gs.lp_sol.size(); i++) {
        if(gs.lp_sol(i)!=0) {
            fmt::print("({}), {}  ", i, gs.lp_sol(i));
        }
    } 
    cout << endl;
    fmt::print("Q2 Gurobi LP status: {}\t Score: {}\n", feasMessage(gs2.checkLpFeasibility(gs2.lp_sol)), gs2.lp_score);
    print_sol(gs2.lp_sol);
    cout << endl << "--------------\n"; 
    fmt::print("Q1 Gurobi ILP status: {}\t Score: {}\n", feasMessage(gs.checkIlpFeasibility(gs.ilp_sol)), gs.ilp_score);
    print_sol(gs.ilp_sol);
    cout << endl;
    fmt::print("Q2 Gurobi ILP status: {}\t Score: {}\n", feasMessage(gs2.checkIlpFeasibility(gs2.ilp_sol)), gs2.ilp_score);
    print_sol(gs2.ilp_sol);
    cout << endl << endl;
    IncrementalEvalGate ieg = IncrementalEvalGate(q1, q2);
    // int cnt=0;
    // VectorXd p = ieg.probe(0);
    // for (int i=0; i<p.size(); i++) {
    //     if(p(i)!=0) {
    //         fmt::print("({}), {}  ", i, p(i));
    //         cnt++;
    //     }
    // }
    // fmt::print("cnt: {}\n", cnt);
    // int feasStatus = gs2.checkIlpFeasibility(p);
    Checker ch = Checker(q2);
    map<long long, long long> sol = ieg.evaluate(10);
    fmt::print("Probe status: {}\t score: {}\nProbing predicted:\n", feasMessage(ch.checkIlpFeasibility(sol)), ch.getScore(sol));
    print_sol(sol);
    sol = ieg.solveSimplifiedILP();
    fmt::print("SimILP status: {}\t score {}\nSimILP predicted:\n", feasMessage(ch.checkIlpFeasibility(sol)), ch.getScore(sol));
    print_sol(sol);

    // cout << feasMessage(gs.checkLpFeasibility(gs2.lp_sol)) << endl;
    // cout << feasMessage(gs2.checkLpFeasibility(gs.lp_sol)) << endl;
    // cout << feasMessage(gs.checkIlpFeasibility(gs2.ilp_sol)) << endl;
    // cout << feasMessage(gs2.checkIlpFeasibility(gs.ilp_sol)) << endl;
    
    return 0;
}

int test_multiple_ILP(){
 vector<int> att_senses = {0, 1, 2}; // Lower bounded, Upper bounded, Both bounded

    // vector<string> cols = {"tmass_prox", "j", "h", "k"};
    // //kt.partitionTable("ssds_3_1", "s500_d5", cols, 500, 5);
    // vector<string> att_cols = {"j", "h", "k"};
    // DetSql det_sql = DetSql("ssds_5_1", "tmass_prox", false, att_cols, att_senses, true, 1);

    vector<string> cols = {"price", "quantity", "discount", "tax"};
    vector<string> att_cols = {"quantity", "discount", "tax"};
    DetSql det_sql = DetSql("tpch_5_1", "price", false, att_cols, att_senses, true, 1);
    
    //det_sql.addFilterWithRatio("j", 0.9, 0);
    //det_sql.addFilterWithRatio("k", 0.9, 2);

    IOFormat printVecFmt(2, 0, ",", " ", "", "", "[", "]");

    // double E = 20; // 50 or 500 for tpc-h
    // double h = 1;
    // LsrProb lsr_prob = LsrProb(det_sql, "s500", 1);

    // q1.display(printVecFmt);
    // q2.display(printVecFmt);

    int num_queries = 10;
    vector<double> h_delta_list = {0.1, 1, 3};
    vector<int> initial_h_list = {1, 3};
    vector<int> E_list = {20, 50};

    DetProb q1 = DetProb(det_sql, -1, 1);
    DetProb q2 = DetProb(det_sql, -1, 1);
    for(auto initial_h : initial_h_list) {
        for(auto h_delta : h_delta_list) {
            for(auto E : E_list) {
                int num_probe = 0;
                int num_similp = 0;
                MeanVar mv_probe, mv_similp;
                for(int i=0; i<num_queries; i++) {
                    q1.generateBounds(E, 1, initial_h);
                    q2.generateBounds(E, 1, initial_h+h_delta);
                    Checker ch = Checker(q2);
                    IncrementalEvalGate ieg = IncrementalEvalGate(q1, q2);
                    map<long long, long long> sol1 = ieg.evaluate(10);
                    map<long long, long long> sol2 = ieg.solveSimplifiedILP();
                    if (ch.checkIlpFeasibility(sol1) == 1) {
                        num_probe++;
                    }
                    if (ch.checkIlpFeasibility(sol2) == 1) num_similp++;
                }
                fmt::print("E: {}\tinitial h: {}\t new h: {}\n", E, initial_h, initial_h+h_delta);
                fmt::print("[Probe]: {}/{}\t[SimILP]: {}/{}\n", num_probe, num_queries, num_similp, num_queries);
            }
        }
    }
    
    return 0;
}

int main() {
    //test_single_ILP();
    test_multiple_ILP();
}