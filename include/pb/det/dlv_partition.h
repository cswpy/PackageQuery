#pragma once

#include "pb/det/lsr_prob.h"

#include "pb/util/upostgres.h"
#include "pb/util/udeclare.h"
#include "pb/util/uconfig.h"
#include "pb/util/udebug.h"
#include "libpq-fe.h"

using namespace pb;

class DLVPartition{
private:
  const LsrProb *prob;
  PgManager *pg;
  PGconn *_conn;
  PGresult *_res;
  string _sql;
public:
  vector<string> cols;
  double group_ratio;
  long long lp_size;
  int layer_count;
  bool is_filtering;
  #if DEBUG
    Profiler pro;
  #endif
private:
  void getNeighboringGroupsRecurse(unordered_set<long long> &group_ids, vector<double> &values, int layer, int index, vector<pair<double, double>> &intervals);
public:
  ~DLVPartition();
  DLVPartition(const LsrProb *prob, vector<string> cols, double group_ratio, long long lp_size, int layer_count);
  string getPName(int layer);
  string getGName(int layer);
  string getInitialGName(int layer);
  bool isCompatible();
  vector<pair<double, double>> getGroupIntervals(int layer, long long group_id);
  long long getGroupContaining(int layer, vector<double> values);
  long long getGroupSize(int layer, long long group_id);
  pair<long long, double> getGroupWorthness(int layer, long long group_id);
  void getGroupComp(vector<long long> &ids, int layer, long long group_id);
  void getNeighboringGroups(unordered_set<long long> &group_ids, int layer, long long group_id);
};