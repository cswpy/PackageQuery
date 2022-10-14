#include "det_exp.h"

#include "pb/util/udeclare.h"
#include "pb/util/uconfig.h"
#include "pb/det/synthetic.h"
#include "pb/det/dlv.h"

using namespace pb;

vector<double> DetExp::H3 = {1, 7, 13};
vector<double> DetExp::H8 = {1, 3, 5, 7, 9, 11, 13, 15};
vector<double> DetExp::E2 = {10, 1000};
vector<double> DetExp::M6 = {8, 16, 32, 64, 128, 300};
vector<double> DetExp::F5 = {0.1, 0.3, 0.5, 0.7, 0.9};

vector<int> DetExp::C6 = {1, 4, 8, 16, 32, 80};
vector<int> DetExp::o4 = {6, 7, 8, 9};

vector<long long> DetExp::N5 = {10000, 100000, 1000000, 10000000, 100000000};

vector<string> DetExp::datasets = {
  "tpch", 
  "ssds"
};

vector<string> DetExp::obj_cols = {
  "price",
  "tmass_prox"
};

vector<bool> DetExp::is_maximizes = {
  true, 
  false
};

vector<vector<string>> DetExp::arr_att_cols = {
  {"quantity", "discount", "tax"},
  {"j", "h", "k"}
};

vector<vector<int>> DetExp::arr_att_senses = {
  {LowerBounded, UpperBounded, Bounded},
  {LowerBounded, UpperBounded, Bounded}
};

vector<bool> DetExp::has_count_constraints = {
  true,
  true
};

vector<long long> DetExp::us = {
  1,
  1
};

DetExp::~DetExp(){
  out.close();
  backup.open(backup_path, std::ios::out);
  for (string s : lines) backup << s;
  backup.close();
  delete pg;
}

DetExp::DetExp(string out_file, bool verbose): verbose(verbose){
  if (verbose) cout << "Start experiment " << out_file << endl;
  string path = fmt::format("{}{}{}{}{}.csv", kProjectHome, separator(), kOutFolder, separator(), out_file);
  backup_path = fmt::format("{}{}{}{}{}{}{}.csv", kProjectHome, separator(), kOutFolder, separator(), kBackupFolder, separator(), out_file);
  out.open(path, std::ios::out);
  pg = new PgManager();
  reset();
}

void DetExp::reset(){
  // Default values
  E = 50;
  a = 0;
  H = 7;
  F = 0.5;
  g = 0.01;
  M = kMainMemorySize;
  S = kLpSize;
  C = kPCore;
  o = 7;
  q = 0;
  seed = 1;
  partition_name = "E0";
}

string DetExp::getTableName(){
  return fmt::format("{}_{}_{}", datasets[q], o, seed);
}

vector<string> DetExp::getCols(){
  vector<string> cols = arr_att_cols[q];
  cols.insert(cols.begin(), obj_cols[q]);
  return cols;
}

DetSql DetExp::generate(){
  string table_name = getTableName();
  if (!pg->existTable(table_name)){
    Synthetic syn = Synthetic();
    syn.createSubtable(datasets[q], o, getCols(), seed);
  }
  DetSql det_sql = DetSql(table_name, obj_cols[q], is_maximizes[q], arr_att_cols[q], arr_att_senses[q], has_count_constraints[q], us[q]);
  return det_sql;
}

void DetExp::partition(){
  DynamicLowVariance dlv = DynamicLowVariance(kPCore, g, kMainMemorySize);
  string table_name = getTableName();
  if (!dlv.existPartition(table_name, partition_name)){
    dlv.partition(table_name, partition_name);
  }
}

void DetExp::write(string id, double x, double y){
  string s = fmt::format("{},{:.2Lf},{:.{}Lf}\n", id, x, y, kPrecision);
  lines.push_back(s);
  if (verbose) cout << s;
  out << s;
  out.flush();
}