#pragma once

#include "pb/util/udebug.h"
#include "pb/util/upostgres.h"
#include "libpq-fe.h"

class DynamicLowVariance{
private:
  string dbname, _sql;
  PgManager *pg;
  PGconn *_conn;
  PGresult *_res;
public:
  static double kGroupRatio;
  static int kMinPartitionSize;
  static int kInitialSize;
  Profiler pro;
private:
  void init();
  void writeStats(string table_name, Stat *stat);
  Stat* readStats(string table_name);
public:
  ~DynamicLowVariance();
  DynamicLowVariance(string dbname);
  void partition(string table_name);
  void partition(string table_name, vector<string> cols);
};