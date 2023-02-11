#pragma once

#include "pb/util/udebug.h"
#include "pb/util/upostgres.h"
#include "libpq-fe.h"

static const double kGroupRatio = 0.01;

class DynamicLowVariance{
private:
  int core;
  double group_ratio, main_memory;
  long long tps;
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
  DynamicLowVariance(int core, double group_ratio=kGroupRatio, double main_memory=kMainMemorySize, long long tps=kLpSize);
  void dropAllPartitions();
  void dropTempTables();
  bool existPartition(string table_name, string partition_name);
  void dropPartition(string table_name, string partition_name);
  void partition(string table_name, string partition_name);
  void partition(string table_name, string partition_name, const vector<string> &cols);
};