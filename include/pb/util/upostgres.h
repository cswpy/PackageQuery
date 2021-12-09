#pragma once

#include "libpq-fe.h"
#include "Eigen/Dense"

using std::string;
using std::vector;
using Eigen::VectorXd;

class Stat{
public:
  long long size;
  VectorXd mean;
  VectorXd M2;
  vector<string> cols;
public:
  Stat(vector<string> cols);
  void add(long long size, VectorXd &mean, VectorXd &M2);
};

class PgManager{
private:
  string dbname, _sql;
  PGconn *_conn;
  PGresult *_res;
public:
  ~PgManager();
  PgManager(string dbname);
  void createExtention(string extension);
  long long getSize(string table_name);
  vector<string> getNumericCols(string table_name);
  Stat* computeStats(string table_name);
  Stat* computeStats(string table_name, const vector<string> &cols);
};