#include <chrono>
#include "pb/util/uconfig.h"
#include "pb/util/upostgres.h"
#include "pb/util/unumeric.h"

#include "kd_tree.h"

KDTree::~KDTree() = default;

KDTree::KDTree(): parent(nullptr), left_child(nullptr), right_child(nullptr) {

}

KDTree::KDTree(int dimension, int size_req, double diameter_req): num_dim(dimension), size_req(size_req), diameter_req(diameter_req) {
    region_mv = MeanVar(dimension);
    parent = nullptr;
    left_child = nullptr;
    right_child = nullptr;
}

// void KDTree::insert(VectorXd tuple) {
//     // Check whether the current quad tree node has capacity
//     // If not, subdivide the node into subgroups and distribute tuples to children nodes
//     tuples.push_back(tuple);
//     mv.add(tuple);
// }

bool KDTree::buildTree(const vector<Tuple>::iterator &begin, const vector<Tuple>::iterator &end, const size_t &length) {
    if (begin==end) {
        return true;
    }

    // cout << distance(begin, end) << " " << length << endl;
    
    // Compute mean/median/var for every subtree
    for (auto it=begin; it!=end; it++) {
        region_mv.add((*it).tup);
    }
    
    // If this is leaf, record all ids
    if((int)length<=size_req && region_mv.getRange().maxCoeff()<=diameter_req) {
        for (auto it=begin; it!=end; it++){
            ids.push_back((*it).id);
        }
        return true;
    }

    // Otherwise, split the tuples into two subtrees
    // Selecting the dimension with largest variance
    
    region_mv.getVar().maxCoeff(&split_axis);

    // Finding the median to split into two subtrees
    auto mid = begin + distance(begin, end)/2;
    nth_element(begin, mid, end, TupleComparator(split_axis));
    median = *mid;
    mid = begin + (length / 2);

    auto left_begin = begin;
    auto right_end = end;

    size_t left_len = length / 2;
    size_t right_len = length - left_len;

    if(left_len > 0) {
        left_child = std::make_shared<KDTree>(num_dim, size_req, diameter_req);
        left_child.get()->buildTree(left_begin, mid, left_len);
    }
    if(right_len > 0) {
        right_child = std::make_shared<KDTree>(num_dim, size_req, diameter_req);
        right_child.get()->buildTree(mid, right_end, right_len);
    }
    return true;
}

void KDTree::getGroupCentroidAndPartition(vector<Tuple>& centroids, map<long long, vector<long long>>& mappings) {
    if(ids.size() > 0){
        Tuple centroid;
        centroid.id = centroids.size();
        centroid.tup = region_mv.getMean();
        centroids.push_back(centroid);
        mappings.insert(pair<long long, vector<long long>>(centroid.id, ids));
    }
    
    if(left_child) left_child.get()->getGroupCentroidAndPartition(centroids, mappings);
    if(right_child) right_child.get()->getGroupCentroidAndPartition(centroids, mappings);
}

void KDTree::partitionTable(string data_table_name, string repr_table_name, vector<string>& cols, int size_req, double diameter_req) {
    PgManager pg = PgManager();
    //long long n_tuples = pg.getSize(table_name);
    int k = cols.size();
    num_dim = k;
    region_mv = MeanVar(k);
    this->size_req = size_req;
    this->diameter_req = diameter_req;
    string col_names = join(cols, ",");
    string query = fmt::format("SELECT {},{} FROM {}", kId, col_names, data_table_name);
    PGconn *conn = PQconnectdb(pg.conninfo.c_str());
    assert(PQstatus(conn) == CONNECTION_OK);
    
    auto t1 = std::chrono::high_resolution_clock::now();
    PGresult *res = NULL;
    res = PQexec(conn, query.c_str());

    vector<Tuple> tuples;
    for(int i=0; i<PQntuples(res); i++) {
        Tuple tuple;
        tuple.id = atol(PQgetvalue(res, i, 0));
        tuple.tup = VectorXd (k);
        for(int j=1; j<=k; j++){
            tuple.tup(j-1) = atof(PQgetvalue(res, i, j));
        }
        tuples.push_back(tuple);
    }
    PQclear(res);
    
    auto t2 = std::chrono::high_resolution_clock::now();
    auto query_elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / 1000000.0;
    fmt::print("Query time elapsed: {:.5Lf}ms\n", query_elapsed);
    buildTree(tuples.begin(), tuples.end(), (size_t)tuples.size());
    auto t3 = std::chrono::high_resolution_clock::now();
    auto build_elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t3 - t2).count() / 1000000.0;
    fmt::print("Finished building KDTree: {:.5Lf}ms\n", build_elapsed);
    
    vector<Tuple> repr_tuples;
    map<long long, vector<long long>> mappings;
    getGroupCentroidAndPartition(repr_tuples, mappings);
    int cnt = 0;
    for(auto repr_tuple : repr_tuples) {
        cnt += mappings[repr_tuple.id].size();
    }
    fmt::print("Created {} groups, with a total of {} tuples\n", repr_tuples.size(), cnt);

    string cid_column_name = "cid";
    if (pg.existTable(repr_table_name)) {
        string alter_table_query = fmt::format(
                                    "ALTER TABLE {0} ADD COLUMN IF NOT EXISTS {1} BIGINT;\n", 
                                    data_table_name, cid_column_name
                                    );

        res = PQexec(conn, alter_table_query.c_str());
        assert(PQresultStatus(res) == PGRES_COMMAND_OK);
        PQclear(res);

        string update_query_template = "UPDATE {0} SET {1}={2}\n"
                                       "WHERE {3} IN ({4});\n";

        // Storing partition info for data table        
        for(auto repr_tuple : repr_tuples) {
            vector<string> sids;
            for(auto id : mappings[repr_tuple.id]) {
                sids.push_back(to_string(id));
            }
            string update_query = fmt::format(
                                    update_query_template, 
                                    data_table_name, 
                                    cid_column_name, 
                                    to_string(repr_tuple.id), 
                                    kId, 
                                    join(sids, ",") 
                                    );
            res = PQexec(conn, update_query.c_str());
            assert(PQresultStatus(res) == PGRES_COMMAND_OK);
            PQclear(res);
        }
        // Truncating the table if there was a stale partition
        string truncate_query = fmt::format("TRUNCATE {}", repr_table_name);
        res = PQexec(conn, truncate_query.c_str());
        assert(PQresultStatus(res) == PGRES_COMMAND_OK);
        PQclear(res);

        // Storing representative info for repr table
        string insert_query_template = "INSERT INTO {0}\n"
                                        "VALUES ({1}, {2});";
        for(auto repr_tuple : repr_tuples) {
            string values = join(repr_tuple.tup, kPrecision);
            string insert_query = fmt::format(insert_query_template, repr_table_name, to_string(repr_tuple.id), values);
            res = PQexec(conn, insert_query.c_str());
            assert(PQresultStatus(res) == PGRES_COMMAND_OK);
            PQclear(res);
        }
    }
    auto t4 = std::chrono::high_resolution_clock::now();
    auto update_elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3).count() / 1000000.0;
    fmt::print("Finished storing partitions: {:.5Lf}ms\n", update_elapsed);

    // Verifying the correctness of partitions
    vector<string> agg_clause, agg_cond;
    for(auto col : cols) {
        agg_clause.push_back(fmt::format("MAX({0})-MIN({0}) AS {0}_range", col));
        agg_cond.push_back(fmt::format("{0}_range>{1}", col, diameter_req));
    }
    string agg_clauses = join(agg_clause, ",\n");

    string test_query_template ="SELECT * FROM (" 
                                "SELECT {0}, COUNT(*) AS size,"
                                "{1}\n"
                                "from {2}\n"
                                "GROUP BY {0} ORDER BY {0}\n"
                                ") p\n"
                                "WHERE {3} OR size>{4}\n";
    string test_query = fmt::format(test_query_template, cid_column_name, agg_clauses, data_table_name, join(agg_cond, " OR "), size_req);
    res = PQexec(conn, test_query.c_str());
    assert(PQntuples(res) == 0);
    fmt::print("Verified paritions satisfy size and diameter requirements\n");
    PQclear(res);
    PQfinish(conn);
}
