#pragma once

#include "pb/util/udeclare.h"
#include "pb/util/unumeric.h"

using namespace pb;

struct Tuple {
    long long id;
    VectorXd tup;
};

struct TupleComparator {
    explicit TupleComparator(int split_axis): split_axis(split_axis) {}

    bool operator () (const Tuple &a, const Tuple &b) {
        return a.tup(split_axis) > b.tup(split_axis);
    }

    int split_axis;
};

class KDTree {
    public:
        int gid, num_dim, size_req, split_axis;
        double diameter_req;
        shared_ptr<KDTree> parent, left_child, right_child;
        Tuple median;
        vector<long long> ids;
        MeanVar region_mv;
    public:
        ~KDTree();
        KDTree();
        KDTree(int dimension, int size_req, double diameter_req);
        void partitionTable(string data_table_name, string partition_name, vector<string>& cols, int size_req, double diameter_req);
        void insert(VectorXd tuple);
        //void storePartitions(string repr_table_name, vector<Tuple>& centroids, map<long long, vector<long long>>& mappings);
        void getGroupCentroidAndPartition(vector<Tuple>& centroids, map<long long, vector<long long>>& mappings, map<long long, double>& diameters);
        bool buildTree(const vector<Tuple>::iterator &begin, const vector<Tuple>::iterator &end, const size_t &length);
};