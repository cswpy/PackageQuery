#pragma once

#include "pb/det/det_prob.h"
#include "pb/det/lsr_prob.h"
#include "pb/det/dlv_partition.h"

#include "pb/util/udebug.h"
#include "pb/util/udeclare.h"
#include "pb/util/upostgres.h"
#include "libpq-fe.h"

using namespace pb;

class LayeredSketchRefine{
private:
  // We need to keep pg as pointer to PgManager since we want to implement customized destructor of LSR and DLV
  PgManager *pg;
  PGconn *_conn;
  PGresult *_res;
  DLVPartition *partition;
  string _sql;
public:
  map<long long, long long> ilp_sol;
  map<long long, double> lp_sol;
  double ilp_score, exe_ilp, lp_score, exe_lp;
  int status;
  Profiler pro;
private:
  void init();
  void formulateDetProb(int core, const LsrProb &prob, DetProb &det_prob, string current_gtable, const vector<long long> &ids);
  DLVPartition* getDLVPartition(const LsrProb *prob);
public:
  ~LayeredSketchRefine();
  LayeredSketchRefine(int core, const LsrProb &prob);
};