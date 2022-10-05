#pragma once

#include "pb/util/udebug.h"
#include "pb/util/upostgres.h"
#include "libpq-fe.h"

class DynamicLowVariance{
private:
  int core;
  double group_ratio;
  string _sql;
  PgManager *pg;
  PGconn *_conn;
  PGresult *_res;
public:
  double exe;
  #if DEBUG
    Profiler pro;
  #endif
private:
  void init();
  long long doPartition(string table_name, string suffix, const vector<string> &cols);
public:
  ~DynamicLowVariance();
  DynamicLowVariance(int core, double group_ratio=0.01);
  void dropAllPartitions();
  void dropTempTables();
  bool existPartition(string table_name, string partition_name);
  void dropPartition(string table_name, string partition_name);
  void partition(string table_name, string partition_name);
  void partition(string table_name, string partition_name, const vector<string> &cols);
};