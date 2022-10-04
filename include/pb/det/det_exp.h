#pragma once

#include "pb/util/upostgres.h"
#include "pb/det/det_sql.h"

class DetExp{
private:
  PgManager *pg;
public:
  // Expected number of tuples in solution package
  double E;
  // Variness of means of constraint bounds
  double a;
  // Hardness of the query
  double H;
  // Filter ratio of the base predicate
  double F;
  // Group ratio of the partition
  double g;
  // Main memory used in partition in GB
  double M;
  // Targeted size of LP in partition
  long long S;
  // Number of cores used
  int C;
  // Order of table
  int o;
  // Query index
  int q;
  // Seed
  int seed;
  // Partition name;
  string partition_name;

  // Ranges
  static vector<double> H3, H8, E2, M6, F5;
  static vector<int> C6, o4;
  static vector<long long> N5;

  // Queries
  static vector<string> datasets;
  static vector<string> obj_cols;
  static vector<bool> is_maximizes;
  static vector<vector<string>> arr_att_cols;
  static vector<vector<int>> arr_att_senses;
  static vector<bool> has_count_constraints;
  static vector<long long> us;
  
public:
  ~DetExp();
  DetExp();
  void reset();
  string getTableName();
  vector<string> getCols();
  DetSql generate();
  void partition();
};
