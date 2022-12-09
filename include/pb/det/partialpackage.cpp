#include "partialpackage.h"
#include "pb/core/gurobi_solver.h"
#include "pb/core/checker.h"

PartialPackage::~PartialPackage() {}

void PartialPackage::init(PGconn *conn, DetProb sketch_det_prob, map<long long, long long> &sketch_sol, unordered_set<long long> &sketch_gids)
{
    this->_conn = conn;
    this->sketch_det_prob = sketch_det_prob;
    this->refine_det_prob = sketch_det_prob;
    this->sketch_sol = sketch_sol;
    this->sketch_gids = sketch_gids;
    this->temp_A = sketch_det_prob.A;
    this->initial_ids = sketch_det_prob.ids;
}

// PartialPackage::PartialPackage(DetProb &det_prob, map<long long, long long> &sketch_sol, vector<int> &group_indices, vector<string> &g_cols, string filter_conds): sketch_sol(sketch_sol), group_indices(group_indices), g_cols(g_cols), filter_conds(filter_conds){
//     num_total_group = det_prob.ids.size();
//     sketch_det_prob = det_prob;
//     refine_det_prob = det_prob;
//     sketch_gids = det_prob.ids;
// }

PartialPackage::PartialPackage(LsrProb &lsr_prob) : lsr_prob(lsr_prob)
{
    g_cols = vector<string>(lsr_prob.det_sql.att_cols.size());
    for (size_t i = 0; i < lsr_prob.det_sql.att_cols.size(); i++)
    {
        g_cols[i] = "g." + lsr_prob.det_sql.att_cols[i];
    }
    filter_conds = getFilterConds(lsr_prob.det_sql.filter_cols, lsr_prob.det_sql.filter_intervals, kPrecision);
}

void PartialPackage::refine(map<long long, long long> &sol)
{
    if (sketch_gids.empty())
        return;

    string g_cols_name = join(g_cols, ",");
    string group_table_name = fmt::format("[1G]_{}_{}", lsr_prob.det_sql.table_name, lsr_prob.partition_name);
    string partition_table_name = fmt::format("[1P]_{}_{}", lsr_prob.det_sql.table_name, lsr_prob.partition_name);

    /* Procedure for refining the groups:
            1. Query Postgres for the actual tuples for a specific group id -> det_prob.A
            2. Remove the group repr tuples from sketch
            3. recalculate the bounds based on other sketch tuples -> det_prob.bl & bu
            4. Copy bounds from lsr_prob for cl, cu, u, l
            5. Solve it using Gurobi

            100 million
            ~10-200 million operations

            1 million = 100 millisecond
    */

    // RMatrixXd effective_A = refine_det_prob.A(Eigen::all, group_indices);
    // auto mean = effective_A.colwise().sum() / effective_A.cols();
    long E = sketch_det_prob.ids.size(); // Package size
    fmt::print("Total package size: {}\n", E);
    int m = g_cols.size() + 1;
    int num_group_refined_phase1 = 0;

    while (sketch_gids.size() > 0)
    {
        // As long as there are unrefined groups, start refining from the first one
        long long refine_group_start_idx = 0, refine_group_end_idx;
        int num_group_refined_iter = 0;
        while (refine_group_start_idx < E)
        {
            long long refine_gid = sketch_det_prob.ids[refine_group_start_idx];
            long long num_repr_tuple = sketch_sol.at(refine_gid);
            refine_group_end_idx = refine_group_start_idx + num_repr_tuple;

            // If this group is already refined, skip
            if (sketch_gids.find(refine_gid) == sketch_gids.end())
            {
                refine_group_start_idx = refine_group_end_idx;
                continue;
            }

            // fmt::print("refine gid: {}\tnum_repr_tuple:{}\trefine_group_start_idx:{}\trefine_group_end_idx:{}\n", refine_gid, num_repr_tuple, refine_group_start_idx, refine_group_end_idx);
            // vector<int> group_indices;
            int *group_idx_seq = new int[E - (long)num_repr_tuple];
            iota(group_idx_seq, group_idx_seq + refine_group_start_idx, 0);
            iota(group_idx_seq + refine_group_start_idx, group_idx_seq + E - num_repr_tuple, refine_group_end_idx);
            // for(int i=0; i<E; i++) {
            //     if(i>=refine_group_start_idx && i<refine_group_end_idx) {
            //         continue;
            //     }
            //     group_indices.push_back(i);
            // }
            // TODO: Handle cl & cu
            RMatrixXd sketch_tmp_A = sketch_det_prob.A(Eigen::all, (*group_idx_seq));
            delete[] group_idx_seq;
            auto sum = sketch_tmp_A.rowwise().sum();
            refine_det_prob.bl = sketch_det_prob.bl - sum;
            refine_det_prob.bu = sketch_det_prob.bu - sum;

            string sql = fmt::format("SELECT {}, {}, {} FROM \"{}\" p INNER JOIN \"{}\" g ON p.tid=g.id WHERE p.gid={} AND {};", kId, lsr_prob.det_sql.obj_col, g_cols_name, partition_table_name, lsr_prob.det_sql.table_name, refine_gid, filter_conds);
            _res = PQexec(_conn, sql.c_str());
            assert(PQresultStatus(_res) == PGRES_TUPLES_OK);

            int n = PQntuples(_res);
            refine_det_prob.resize(m, n);

            for (int i = 0; i < n; i++)
            {
                refine_det_prob.ids[i] = atoll(PQgetvalue(_res, i, 0));
                refine_det_prob.c[i] = atof(PQgetvalue(_res, i, 1));
                for (int j = 0; j < m - 1; j++)
                {
                    refine_det_prob.A(j, i) = atof(PQgetvalue(_res, i, 2 + j));
                }
                if (m > 0)
                    refine_det_prob.A(m - 1, i) = 1.0; // Only change for cl & cu when there is bl & bu
            }
            PQclear(_res);

            refine_det_prob.u.fill((double)lsr_prob.det_sql.u);
            refine_det_prob.truncate();
            n = refine_det_prob.A.cols();
            m = refine_det_prob.A.rows();

            // Solve refine for each group
            GurobiSolver gs = GurobiSolver(refine_det_prob, true);
            gs.solveIlp();
            // cout << "Refine status:" << gs.ilp_status << ": " << solMessage(gs.ilp_status) << endl;

            Checker ch = Checker(refine_det_prob);
            int feasStatus = ch.checkIlpFeasibility(gs.ilp_sol);
            if (feasStatus != 1)
            {
                fmt::print("Refine failed: {} at group {}", feasMessage(feasStatus), refine_gid);
                continue;
            }
            else
            {
                num_group_refined_iter++;
                assert(sketch_gids.erase(refine_gid) > 0);
            }

            // Replace the repr tuples with the actual tuples
            int tmp = refine_group_start_idx;
            for (int i = 0; i < gs.ilp_sol.size() && tmp < refine_group_end_idx; i++)
            {
                if (gs.ilp_sol(i) != 0)
                {
                    // int occurrence = (int)gs.ilp_sol(i);
                    for (int j = 0; j < m; j++)
                    {
                        this->temp_A(j, tmp) = refine_det_prob.A(j, i);
                        // sketch_det_prob.A(j, tmp) = refine_det_prob.A(j, i);
                    }
                    sketch_det_prob.ids[tmp] = refine_det_prob.ids[i]; // Assigning the actual tuples' IDs to sketch_det_prob does not affect the calculation of ILP constraints
                    tmp += 1;
                }
            }
            num_group_refined_phase1 += 1;
            refine_group_start_idx = refine_group_end_idx;
        }

        if (num_group_refined_iter == 0)
        {
            fmt::print("[Failed] No group could be refined in this iteration\n");
            return;
        }

        // Update the sketch_det_prob to actual tuples refined in this round
        sketch_det_prob.A = temp_A;
    }

    for (auto id : sketch_det_prob.ids)
    {
        auto pair = sol.emplace(id, 0);
        pair.first->second += 1;
    }

    LsrChecker lsr_ch = LsrChecker(lsr_prob);
    int feasStatus = lsr_ch.checkIlpFeasibility(sol);
    if (feasStatus == 1)
    {
        fmt::print("[Suceed] Phase 1 feasible, finished refining for {} groups, with {} tuples in the package\n", num_group_refined_phase1, sketch_det_prob.ids.size());
        return;
    }

    fmt::print("Phase 1 failed, feasibility status: {}, trying to re-refine\n", feasMessage(feasStatus));

    // Phase 2: Going though each group again for replacement
    long long refine_group_start_idx = 0, refine_group_end_idx;
    while (refine_group_start_idx < E)
    {
        // long long refine_gid = sketch_det_prob.ids[refine_group_start_idx];
        long long re_refine_gid = initial_ids[refine_group_start_idx];
        long long num_repr_tuple = sketch_sol[re_refine_gid];
        refine_group_end_idx = refine_group_start_idx + num_repr_tuple;

        // fmt::print("refine gid: {}\tnum_repr_tuple:{}\trefine_group_start_idx:{}\trefine_group_end_idx:{}\n", refine_gid, num_repr_tuple, refine_group_start_idx, refine_group_end_idx);
        // vector<int> group_indices;
        int *group_idx_seq = new int[E - (long)num_repr_tuple];
        iota(group_idx_seq, group_idx_seq + refine_group_start_idx, 0);
        iota(group_idx_seq + refine_group_start_idx, group_idx_seq + E - num_repr_tuple, refine_group_end_idx);
        // for(int i=0; i<E; i++) {
        //     if(i>=refine_group_start_idx && i<refine_group_end_idx) {
        //         continue;
        //     }
        //     group_indices.push_back(i);
        // }
        // TODO: Handle cl & cu
        RMatrixXd sketch_tmp_A = sketch_det_prob.A(Eigen::all, (*group_idx_seq));
        delete[] group_idx_seq;
        auto sum = sketch_tmp_A.rowwise().sum();
        refine_det_prob.bl = sketch_det_prob.bl - sum;
        refine_det_prob.bu = sketch_det_prob.bu - sum;

<<<<<<< HEAD
        // TODO: check sql syntax
        string sql = fmt::format("SELECT {}, {}, {} FROM \"{}\" p INNER JOIN \"{}\" g ON p.tid=g.id WHERE p.gid={} AND {};", kId, lsr_prob.det_sql.obj_col, g_cols_name, partition_table_name, lsr_prob.det_sql.table_name, refine_gid, filter_conds);
=======
        string sql = fmt::format("SELECT {}, {}, {} FROM \"{}\" p INNER JOIN \"{}\" g ON p.tid=g.id WHERE p.gid={} AND {};", kId, lsr_prob.det_sql.obj_col, g_cols_name, partition_table_name, lsr_prob.det_sql.table_name, re_refine_gid, filter_conds);
>>>>>>> Implemented iterative refine
        _res = PQexec(_conn, sql.c_str());
        assert(PQresultStatus(_res) == PGRES_TUPLES_OK);

        int n = PQntuples(_res);
        refine_det_prob.resize(m, n);

        for (int i = 0; i < n; i++)
        {
            refine_det_prob.ids[i] = atoll(PQgetvalue(_res, i, 0));
            refine_det_prob.c[i] = atof(PQgetvalue(_res, i, 1));
            for (int j = 0; j < m - 1; j++)
            {
                refine_det_prob.A(j, i) = atof(PQgetvalue(_res, i, 2 + j));
            }
            if (m > 0)
                refine_det_prob.A(m - 1, i) = 1.0; // Only change for cl & cu when there is bl & bu
        }
        PQclear(_res);

        refine_det_prob.u.fill((double)lsr_prob.det_sql.u);
        refine_det_prob.truncate();
        n = refine_det_prob.A.cols();
        m = refine_det_prob.A.rows();

        // Solve refine for each group
        GurobiSolver gs = GurobiSolver(refine_det_prob, true);
        gs.solveIlp();
        // cout << "Refine status:" << gs.ilp_status << ": " << solMessage(gs.ilp_status) << endl;

        Checker ch = Checker(refine_det_prob);
        int feasStatus = ch.checkIlpFeasibility(gs.ilp_sol);
        if (feasStatus != 1)
        {
            //fmt::print("Refine failed {} at group {}", feasMessage(feasStatus), refine_gid);
            refine_group_start_idx = refine_group_end_idx;
            continue;
        }

        // Replace the repr tuples with the actual tuples
        int tmp = refine_group_start_idx;
        for (int i = 0; i < gs.ilp_sol.size() && tmp < refine_group_end_idx; i++)
        {
            if (gs.ilp_sol(i) != 0)
            {
                // int occurrence = (int)gs.ilp_sol(i);
                for (int j = 0; j < m; j++)
                {
                    sketch_det_prob.A(j, tmp) = refine_det_prob.A(j, i);
                }
                sketch_det_prob.ids[tmp] = refine_det_prob.ids[i];
                sol.erase(sketch_det_prob.ids[tmp]);
                auto pair = sol.emplace(refine_det_prob.ids[i], 0);
                pair.first->second += 1;
                tmp += 1;
            }
        }
        // Solution updated, return
        fmt::print("[Suceed] Phase 2 suceeded, re-refined group {}, feasMessage: {}\n", re_refine_gid, feasMessage(feasStatus));
        return;

    }

    // Exited the loop means that no re-refining can be done on any group
    fmt::print("[Failed] Phase 2 failed, no re-refining can be done.\n");
    return;
}
