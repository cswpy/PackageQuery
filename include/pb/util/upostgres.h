#pragma once

#include "pb/util/udeclare.h"
#include "pb/util/uconfig.h"

#include "libpq-fe.h"
#include "Eigen/Dense"

using namespace pb;

using std::string;
using std::vector;
using Eigen::VectorXd;

void showError(PGconn *conn);

class Stat{
public:
  long long size;
  VectorXd mean, M2, amin, amax;
  vector<string> cols;
public:
  Stat(vector<string> cols);
  void add(long long size, VectorXd mean, VectorXd M2, VectorXd amin, VectorXd amax);
  void add(long long size, VectorXd mean, VectorXd M2);
  int getIndex(string col);
  double getVar(int i);
  VectorXd getVars();
};

class PgManager{
private:
  string _sql;
  PGconn *_conn;
  PGresult *_res;
private:
  bool checkStats(string table_name);
  void writeStats(string table_name, Stat *stat);
  Stat* computeStats(string table_name, const vector<string> &cols);
public:
  string conninfo;
  ~PgManager();
  PgManager();
  long long getSize(string table_name);
  vector<string> listTables(string schema_name=kPgSchema);
  vector<string> listColumns(string table_name);
  bool existTable(string table_name);
  void dropTable(string table_name);
  vector<string> getNumericCols(string table_name);
  Stat* readStats(string table_name);
  void deleteStats(string table_name);
  void getSelectedTuples(RMatrixXd &out_tuples, vector<long long> &out_ids, string table_name, vector<long long> &ids, vector<string> cols, vector<string> filter_cols=vector<string>(), vector<pair<double, double>> filter_intervals=vector<pair<double, double>>());
  void getConsecutiveTuples(RMatrixXd &out_tuples, vector<long long> &out_ids, string table_name, long long start_id, long long end_id, vector<string> cols, vector<string> filter_cols=vector<string>(), vector<pair<double, double>> filter_intervals=vector<pair<double, double>>());
};