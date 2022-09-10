#pragma once

#include "libpq-fe.h"
#include "Eigen/Dense"

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
};

class PgManager{
private:
  string _sql;
  PGconn *_conn;
  PGresult *_res;
public:
  string conninfo;
  ~PgManager();
  PgManager();
  long long getSize(string table_name);
  bool existTable(string table_name);
  void dropTable(string table_name);
  vector<string> getNumericCols(string table_name);
  Stat* computeStats(string table_name, const vector<string> &cols);
};