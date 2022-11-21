#include "sr.h"

#include "pb/core/gurobi_solver.h"
#include "pb/det/det_prob.h"
#include "partialpackage.h"

#include "pb/util/uconfig.h"
#include "pb/util/upostgres.h"
#include "pb/util/unumeric.h"

SketchRefine::~SketchRefine()
{
    PQfinish(_conn);
    delete pg;
}

void SketchRefine::init()
{
    pg = new PgManager();
    _conn = PQconnectdb(pg->conninfo.c_str());
    assert(PQstatus(_conn) == CONNECTION_OK);
    _res = NULL;
}

SketchRefine::SketchRefine(LsrProb &lsr_prob): prob(lsr_prob)
{
    init();
}

void SketchRefine::refine(map<long long, long long> &sol) {
    group_table_name = fmt::format("[1G]_{}_{}", prob.det_sql.table_name, prob.partition_name);
    partition_table_name = fmt::format("[1P]_{}_{}", prob.det_sql.table_name, prob.partition_name);
    assert(pg->existTable(group_table_name));
    assert(pg->existTable(partition_table_name));

    // Filtering groups based on base predicates
    if (prob.det_sql.isFiltering())
    {
        long long num_group = pg->getSize(group_table_name);
        vector<string> g_cols(prob.det_sql.att_cols.size());
        for (size_t i = 0; i < prob.det_sql.att_cols.size(); i++){
            g_cols[i] = "g." + prob.det_sql.att_cols[i];
        }
        string g_cols_name = join(g_cols, ",");
        string filter_conds = getFilterConds(prob.det_sql.filter_cols, prob.det_sql.filter_intervals, kPrecision);
        // Construct lower & upper bound for filtering
        // assert(g_cols.size() == prob.det_sql.filter_cols.size());
        VectorXd bl(g_cols.size());
        bl.fill(-DBL_MAX);
        VectorXd bu(g_cols.size());
        bu.fill(DBL_MAX);
        for (size_t i = 0; i < prob.det_sql.filter_cols.size(); i++)
        {
            auto itr = find(prob.det_sql.att_cols.begin(), prob.det_sql.att_cols.end(), prob.det_sql.filter_cols[i]);
            assert(itr != prob.det_sql.att_cols.end());
            int ind = itr - prob.det_sql.att_cols.begin();
            bl(ind) = prob.det_sql.filter_intervals[i].first;
            bu(ind) = prob.det_sql.filter_intervals[i].second;
        }
        // Updating group averages according to filtered tuples
        vector<vector<double>> group_averages;
        group_averages.reserve(num_group);
        vector<long long> group_ids;
        //group_average_flattened.reserve(num_group * g_cols.size());
        for (long long i = 1; i <= num_group; i++)
        {
            string sql = fmt::format("SELECT {} FROM \"{}\" p INNER JOIN \"{}\" g ON p.tid=g.id WHERE p.gid={}{};", g_cols_name, partition_table_name, prob.det_sql.table_name, i, filter_conds);
            _res = PQexec(_conn, sql.c_str());
            assert(PQresultStatus(_res) == PGRES_TUPLES_OK);
            MeanVar group_mv (g_cols.size());
            VectorXd tuple(g_cols.size());
            for (int j = 0; j < PQntuples(_res); j++)
            {
                for (size_t k = 0; k < g_cols.size(); k++) 
                {
                    tuple(k) = atof(PQgetvalue(_res, j, k));
                }
                if (checkTupleFiltered(tuple, bl, bu))
                {
                    group_mv.add(tuple);
                }
            }
            PQclear(_res);
            // add group to Sketch
            if (group_mv.sample_count > 0)
            {
                group_ids.push_back(i);
                VectorXd mean = group_mv.getMean();
                vector<double> group_average (mean.begin(), mean.end());
                group_averages.push_back(group_average);
            }
        }
        int m = g_cols.size() + 1; // Account for cl & cu
        int n = group_ids.size();

        det_prob.resize(m, n);
        for(int i=0; i<n; i++) {
            for (int j=0; j<m-1; j++) {
                det_prob.A(j, i) = group_averages[i][j];
            }
            det_prob.A(m-1, i) = 1.0; // for cl & cu
        }
        for(int i=0; i<n; i++) {
            det_prob.ids[i] = group_ids[i];
        }
    }
    // If no filtering, we use the saved partition data directly
    else
    {
        string col_names = join(prob.det_sql.att_cols, ",");
        string sql = fmt::format("SELECT {},{},{} FROM \"{}\";",
                                 kId, prob.det_sql.obj_col, col_names, group_table_name);
        _res = PQexec(_conn, sql.c_str());
        assert(PQresultStatus(_res) == PGRES_TUPLES_OK);
        int m = prob.det_sql.att_cols.size() + 1;
        int n = PQntuples(_res);
        det_prob.resize(m, n);
        for (int i = 0; i < PQntuples(_res); i++)
        {
            det_prob.ids[i] = atoll(PQgetvalue(_res, i, 0));
            det_prob.c(i) = atof(PQgetvalue(_res, i, 1));
            for (int j = 0; j < m - 1; j++)
            {
                det_prob.A(j, i) = atof(PQgetvalue(_res, i, 2 + j));
            }
            det_prob.A(m - 1, i) = 1.0; // for cl & cu
        }
        PQclear(_res);
    }

    // Assign other common variables
    det_prob.u.fill((double)prob.det_sql.u); // Each tuple's occurrence is bounded between 0 and u
    det_prob.copyBounds(prob.bl, prob.bu, prob.cl, prob.cu);
    det_prob.truncate();

    // Sketching initial package from representative tuples
    GurobiSolver gs = GurobiSolver(det_prob, true);
    gs.solveIlp();

    // Replacing representative tuples with actual tuples
    map<long long, long long> sketch_sol;
    vector<int> sketch_gids; // Different from tid, the indices within a group
    vector<long long> temp_ids;

    //group_indices.reserve((int)prob.cu); 
    for (int i = 0; i < gs.ilp_sol.size(); i++)
    {
        if (gs.ilp_sol(i) != 0)
        {
            sketch_sol.insert(pair<long long, long long>(det_prob.ids[i], gs.ilp_sol(i)));
            //group_indices.push_back(i);
            sketch_gids.insert(sketch_gids.end(), (int)gs.ilp_sol(i), i);
            temp_ids.insert(temp_ids.end(), (int)gs.ilp_sol(i), det_prob.ids[i]);
        }
    }
    RMatrixXd effective_A = det_prob.A(Eigen::all, sketch_gids);
    VectorXd effective_c = det_prob.c(sketch_gids);
    Checker ch = Checker(det_prob);
    fmt::print("Sketch solution status: {}\n",feasMessage(ch.checkIlpFeasibility(sketch_sol)) );

    // Effectively, only changes variables whose dimension includes n
    int m = effective_A.rows();
    int n = effective_A.cols();
    det_prob.resize(m, n);
    det_prob.A = effective_A;
    det_prob.ids = temp_ids;
    det_prob.c = effective_c;

    PartialPackage pp = PartialPackage(prob);
    pp.init(_conn, det_prob, sketch_sol, sketch_gids);
    pp.refine(sol);
}

bool SketchRefine::checkTupleFiltered(VectorXd &tuple, VectorXd &bl, VectorXd &bu)
{   
    return (tuple.array() > bl.array()).all() && (tuple.array() < bu.array()).all();
}

// void SketchRefine::Sketch(LsrProb &prob)
// {
//     string col_names = join(prob.det_sql.att_cols, ",");
//     string sql = fmt::format("SELECT {},{},{} FROM \"{}\";",
//                              kId, prob.det_sql.obj_col, col_names, group_table_name);
//     formulateDetProb(prob, sql);
// }

// void SketchRefine::Refine(LsrProb &prob, long refine_gid)
// {
//     string col_names = join(prob.det_sql.att_cols, ",");
//     string group_select_template = "SELECT {},{},{} FROM \"{}\" d\n"
//                                    "JOIN {} p\n ON d.{}=p.tid\n"
//                                    "WHERE gid={}";
//     string sql = fmt::format(group_select_template,
//                              kId, prob.det_sql.obj_col, col_names, prob.det_sql.table_name,
//                              partition_table_name, kId,
//                              refine_gid);
//     formulateDetProb(prob, sql);
// }

// void SketchRefine::formulateDetProb(LsrProb &prob, string sql)
// {
//     // Execute the filtering query to get base relations
//     _res = PQexec(_conn, sql.c_str());

//     int m = prob.det_sql.att_cols.size() + 1; // number of attributes + 1 for cl & cu
//     int n = PQntuples(_res);                  // number of variables/tuples
//     det_prob.resize(m, n);
//     det_prob.u.fill((double)prob.det_sql.u); // Each tuple's occurrence is bounded between 0 and u
//     det_prob.copyBounds(prob.bl, prob.bu, prob.cl, prob.cu);

//     for (int i = 0; i < PQntuples(_res); i++)
//     {
//         det_prob.ids[i] = atoll(PQgetvalue(_res, i, 0));
//         det_prob.c(i) = atof(PQgetvalue(_res, i, 1));
//         for (int j = 0; j < m - 1; j++)
//         {
//             det_prob.A(j, i) = atof(PQgetvalue(_res, i, 2 + j));
//         }
//         det_prob.A(m - 1, i) = 1.0; // for cl & cu
//     }
//     PQclear(_res);
//     det_prob.truncate();

//     cout << "----bl----" << endl;
//     cout << det_prob.bl << endl;
//     cout << "----bu----" << endl;
//     cout << det_prob.bu << endl;
//     // cout << "----c(obj_col)----" << endl;
//     // cout << det_prob.c << endl;
//     // cout << "----A(data_col)----"<< endl;
//     // cout << det_prob.A << endl;
//     // cout << "----l----" << endl;
//     // cout << det_prob.l << endl;
//     // cout << "----u----" << endl;
//     // cout << det_prob.u << endl;
// }