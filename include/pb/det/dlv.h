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
  static double kGroupRatio, kVarScale;
  static int kMaxSize, kInitialSize;
  Profiler pro;
private:
  void init();
  void writeStats(string table_name, Stat *stat);
  Stat* readStats(string table_name);
  void quickPartition(string table_name, Stat *stat, const vector<string> &cols);
  void doPartition(string table_name, string suffix, const vector<string> &cols);
public:
  ~DynamicLowVariance();
  DynamicLowVariance(string dbname);
  void partition(string table_name, string partition_name);
  void partition(string table_name, string partition_name, const vector<string> &cols);
};