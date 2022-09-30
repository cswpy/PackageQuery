#pragma once

#include "pb/det/det_prob.h"
#include "pb/det/lsr_prob.h"

#include "pb/util/udebug.h"
#include "pb/util/udeclare.h"
#include "pb/util/upostgres.h"
#include "libpq-fe.h"

using namespace pb;

class DLVPartition{
private:
  PgManager *pg;
  PGconn *_conn;
  PGresult *_res;
  string _sql;
public:
  string table_name, partition_name;
  vector<string> cols;
  double group_ratio;
  long long lp_size;
  int layer_count;
public:
  ~DLVPartition();
  DLVPartition(string table_name, string partition_name, vector<string> cols, double group_ratio, long long lp_size, int layer_count);
  string getPName(int layer);
  string getGName(int layer, bool is_filtered=false);
  bool isCompatible(const LsrProb &prob);
  long long getGroupSize(int layer, long long id, bool is_filtered=false);
};

class LayeredSketchRefine{
private:
  // We need to keep pg as pointer to PgManager since we want to implement customized destructor of LSR and DLV
  PgManager *pg;
  PGconn *_conn;
  PGresult *_res;
  DLVPartition *partition;
  string _sql;
public:
  unordered_map<long long, long long> ilp_sol;
  unordered_map<long long, double> lp_sol;
  double ilp_score, exe_ilp, lp_score, exe_lp;
  int status;
  Profiler pro;
private:
  void init();
  void formulateDetProb(int core, const LsrProb &prob, DetProb &det_prob, string current_gtable, const vector<long long> &ids);
  DLVPartition* getDLVPartition(string table_name, string partition_name);
public:
  ~LayeredSketchRefine();
  LayeredSketchRefine(int core, const LsrProb &prob);
};