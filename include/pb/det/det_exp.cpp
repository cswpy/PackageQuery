#include "det_exp.h"

#include "pb/util/udeclare.h"
#include "pb/util/uconfig.h"
#include "pb/det/synthetic.h"
#include "pb/det/dlv.h"
#include "pb/det/kd_tree.h"

using namespace pb;

const double kTauRatio = 0.001;

vector<double> DetExp::H8 = {1, 3, 5, 7, 9, 11, 13, 15};
vector<double> DetExp::H4 = {1, 3, 5, 7};
vector<double> DetExp::E2 = {20, 1000};
vector<double> DetExp::M6 = {8, 16, 32, 64, 128, 300};
vector<double> DetExp::F5 = {0.1, 0.3, 0.5, 0.7, 0.9};
vector<double> DetExp::g2 = {0.001, 0.01};

vector<int> DetExp::C7 = {1, 2, 4, 8, 16, 32, 80};
vector<int> DetExp::o6 = {4, 5, 6, 7, 8, 9};
vector<int> DetExp::o4 = {6, 7, 8, 9};

vector<long long> DetExp::S3 = {100000, 1000000, 10000000};

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

vector<string> DetExp::filtered_cols = {
  "price",
  "tmass_prox"
};

vector<double> DetExp::Es = {
  30.0, 30.0
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
  M = kMainMemorySize;
  /****/
  g = kGroupRatio;
  S = kLpSize;
  /****/
  C = kPCore;
  o = 8;
  q = 0;
  tps = kTps;
  tau_ratio = kTauRatio;
  seed = 1;
}

string DetExp::getTableName(){
  return getSubtableName(datasets[q], o, seed);
}

vector<string> DetExp::getCols(){
  vector<string> cols = arr_att_cols[q];
  cols.insert(cols.begin(), obj_cols[q]);
  return cols;
}

DetSql DetExp::generate(bool is_lazy){
  string table_name = getTableName();
  Synthetic syn = Synthetic();
  if (is_lazy){
    if (!pg->existTable(table_name)){
      syn.createSubtable(datasets[q], o, getCols(), seed);
    }
  } else{
    pg->dropTable(table_name);
    syn.createSubtable(datasets[q], o, getCols(), seed);
  }

  DetSql det_sql = DetSql(table_name, obj_cols[q], is_maximizes[q], arr_att_cols[q], arr_att_senses[q], has_count_constraints[q], us[q]);
  return det_sql;
}

double DetExp::dlvPartition(bool is_lazy){
  DynamicLowVariance dlv = DynamicLowVariance(C, g, M, tps);
  string table_name = getTableName();
  if (is_lazy){
    if (!dlv.existPartition(table_name, getDlvPartitionName())){
      dlv.dropPartition(table_name, getDlvPartitionName());
      dlv.partition(table_name, getDlvPartitionName());
    }
  } else{
    dlv.dropPartition(table_name, getDlvPartitionName());
    dlv.partition(table_name, getDlvPartitionName());
  }
  return dlv.exe;
}

double DetExp::kdPartition(bool is_lazy){
  KDTree kt;
  string ptable = fmt::format("[1P]_{}_{}", getTableName(), getKdPartitionName());
  string gtable = fmt::format("[1G]_{}_{}", getTableName(), getKdPartitionName());
  double tau = kTauRatio * pg->getSize(getTableName());
  if (is_lazy){
    if (!pg->existTable(ptable) || !pg->existTable(gtable)){
      pg->dropTable(ptable);
      pg->dropTable(gtable);
      kt.partitionTable(getTableName(), getKdPartitionName(), getCols(), tau, DBL_MAX);
    }
  } else{
    pg->dropTable(ptable);
    pg->dropTable(gtable);
    kt.partitionTable(getTableName(), getKdPartitionName(), getCols(), tau, DBL_MAX);
  }
  return kt.exec_kd;
}

string DetExp::getDlvPartitionName(){
  return fmt::format("dlv_C{}_g{}_M{}_tps{}", C, formatFloat(g), formatFloat(M, 1), tps);
}

string DetExp::getKdPartitionName(){
  return fmt::format("kd_tau{}", formatFloat(tau_ratio));
}

void DetExp::write(string id, double x, double y){
  string s = fmt::format("{},{:.{}Lf},{:.{}Lf}\n", id, x, kPrecision, y, kPrecision);
  lines.push_back(s);
  if (verbose) cout << s;
  out << s;
  out.flush();
}

void DetExp::write(string label, double v){
  string s = fmt::format("{},{:.{}Lf}\n", label, v, kPrecision);
  lines.push_back(s);
  if (verbose) cout << s;
  out << s;
  out.flush();
}