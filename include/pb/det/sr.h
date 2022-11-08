#pragma once

#include "pb/det/det_prob.h"
#include "pb/det/lsr_prob.h"
#include "pb/det/kd_tree.h"

#include "pb/util/udebug.h"
#include "pb/util/udeclare.h"
#include "pb/util/upostgres.h"
#include "pb/util/uconfig.h"
#include "libpq-fe.h"

using namespace pb;

class SketchRefine {
    private:
        PgManager *pg;
        PGconn *_conn;
        PGresult *_res;
        string _sql, group_table_name, partition_table_name;
        DetProb det_prob;
    public:
        ~SketchRefine();
        SketchRefine(LsrProb &lsr_prob);
        void init();
        void Sketch(LsrProb &prob);
        void Refine(LsrProb &prob, long refine_gid);
        void formulateDetProb(LsrProb &prob, string sql);
        bool checkTupleFiltered(VectorXd &tuple, VectorXd &bl, VectorXd &bu);
};
