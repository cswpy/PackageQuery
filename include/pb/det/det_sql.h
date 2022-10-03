#pragma once

#include "pb/util/udeclare.h"
#include "pb/util/upostgres.h"

using namespace pb;

static const double kLowerCountFactor = 0.5;
static const double kUpperCountFactor = 1.5;

class DetSql{
private:
  PgManager *pg;
  vector<string> table_cols;
public:
  string table_name, obj_col;
  bool is_maximize;
  vector<string> att_cols;
  vector<int> att_senses;
  // Assuming upper bounded by u and lower bounded by 0 for all tuples.
  long long u;
  bool has_count_constraint;

  vector<string> filter_cols;
  vector<pair<double, double>> filter_intervals;
public:
  ~DetSql();
  DetSql(string table_name, string obj_col, bool is_maximize, vector<string> att_cols, vector<int> att_senses, bool has_count_constraint=true, long long u=1);
  void addFilter(string col, double l, double u);
  bool isFiltering();
};