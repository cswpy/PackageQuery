#include "dlv_partition.h"

#include "pb/util/uconfig.h"

const double kIntervalEps = 1e-6;

DLVPartition::~DLVPartition(){
  delete pg;
  _sql = "DEALLOCATE ALL;";
  _res = PQexec(_conn, _sql.c_str());
  PQclear(_res);
  PQfinish(_conn);
}

DLVPartition::DLVPartition(const LsrProb *prob, vector<string> cols, double group_ratio, long long lp_size, int layer_count)
  : prob(prob), cols(cols), group_ratio(group_ratio), lp_size(lp_size), layer_count(layer_count){
  pg = new PgManager();
  _conn = PQconnectdb(pg->conninfo.c_str());
  assert(PQstatus(_conn) == CONNECTION_OK);
  _res = NULL;

  is_filtering = prob->det_sql.filter_cols.size() > 0;

  vector<string> intervals (cols.size());
  for (int i = 0; i < (int) cols.size(); i ++) intervals[i] = "interval_" + cols[i];
  string interval_names = join(intervals, ",");

  vector<string> itv_conds (cols.size());
  for (int i = 0; i < (int) cols.size(); i ++) itv_conds[i] = fmt::format("${}::float8 <@ {}", i+1, intervals[i]);
  string and_itv_conds = join(itv_conds, " AND ");

  for (int layer=1; layer <= layer_count; layer ++){
    string current_gtable = getGName(layer);
    _sql = fmt::format("SELECT size FROM \"{}\" WHERE id=$1::bigint;", current_gtable);
    _res = PQprepare(_conn, fmt::format("size_{}", current_gtable).c_str(), _sql.c_str(), 1, NULL);
    assert(PQresultStatus(_res) == PGRES_COMMAND_OK);
    PQclear(_res);
    _sql = fmt::format("SELECT size,{} FROM \"{}\" WHERE id=$1::bigint;", prob->det_sql.obj_col, current_gtable);
    _res = PQprepare(_conn, fmt::format("worth_{}", current_gtable).c_str(), _sql.c_str(), 1, NULL);
    assert(PQresultStatus(_res) == PGRES_COMMAND_OK);
    PQclear(_res);

    string current_initial_gtable = getInitialGName(layer);
    _sql = fmt::format("SELECT {} FROM \"{}\" WHERE id=$1::bigint;", interval_names, current_initial_gtable);
    _res = PQprepare(_conn, fmt::format("interval_{}", current_initial_gtable).c_str(), _sql.c_str(), 1, NULL);
    assert(PQresultStatus(_res) == PGRES_COMMAND_OK);
    PQclear(_res);
    _sql = fmt::format("SELECT id FROM \"{}\" WHERE {};", current_initial_gtable, and_itv_conds);
    _res = PQprepare(_conn, fmt::format("group_{}", current_initial_gtable).c_str(), _sql.c_str(), cols.size(), NULL);
    assert(PQresultStatus(_res) == PGRES_COMMAND_OK);
    PQclear(_res);

    string current_ptable = getPName(layer);
    if (!is_filtering){
      _sql = fmt::format("SELECT tid FROM \"{}\" WHERE gid=$1::bigint;", current_ptable);
    } else{
      if (layer == 1){
        // Special filtering for original layer
        string filter_conds = getFilterConds(prob->det_sql.filter_cols, prob->det_sql.filter_intervals, kPrecision);
        _sql = fmt::format("SELECT p.tid FROM \"{}\" p INNER JOIN \"{}\" g ON p.tid=g.id WHERE p.gid=$1::bigint{};", current_ptable, prob->det_sql.table_name, filter_conds);
      } else{
        _sql = fmt::format("SELECT p.tid FROM \"{}\" p WHERE p.gid=$1::bigint AND p.tid IN (SELECT id FROM \"{}\");", current_ptable, getGName(layer-1));
      }
    }
    _res = PQprepare(_conn, fmt::format("comp_{}", current_ptable).c_str(), _sql.c_str(), 1, NULL);
    assert(PQresultStatus(_res) == PGRES_COMMAND_OK);
    PQclear(_res);
  }
}

string DLVPartition::getPName(int layer){
  return fmt::format("[{}P]{}_{}", layer, prob->det_sql.table_name, prob->partition_name);
}

string DLVPartition::getGName(int layer){
  if (!is_filtering) return fmt::format("[{}G]{}_{}", layer, prob->det_sql.table_name, prob->partition_name);
  else return fmt::format("{}_[{}G]{}_{}", kTempPrefix, layer, prob->det_sql.table_name, prob->partition_name);
}

string DLVPartition::getInitialGName(int layer){
  return fmt::format("[{}G]{}_{}", layer, prob->det_sql.table_name, prob->partition_name);
}

bool DLVPartition::isCompatible(){
  bool res = isIn(cols, prob->det_sql.obj_col);
  for (string s : prob->det_sql.att_cols) res &= isIn(cols, s);
  return res;
}

long long DLVPartition::getGroupSize(int layer, long long group_id){
  string current_gtable = getGName(layer);
  char* vals[1];
  assign(vals, 0, group_id);
  _res = PQexecPrepared(_conn, fmt::format("size_{}", current_gtable).c_str(), 1, vals, NULL, NULL, 0);
  free(vals, 1);
  long long size = 0;
  if (PQntuples(_res)){
    size = atoll(PQgetvalue(_res, 0, 0));
  }
  PQclear(_res);
  return size;
}

pair<long long, double> DLVPartition::getGroupWorthness(int layer, long long group_id){
  string current_gtable = getGName(layer);
  char* vals[1];
  assign(vals, 0, group_id);
  _res = PQexecPrepared(_conn, fmt::format("worth_{}", current_gtable).c_str(), 1, vals, NULL, NULL, 0);
  free(vals, 1);
  long long size = 0;
  double worth = -1;
  if (PQntuples(_res)){
    size = atoll(PQgetvalue(_res, 0, 0));
    worth = atof(PQgetvalue(_res, 0, 1));
    if (!prob->det_sql.is_maximize) worth = -worth;
  }
  PQclear(_res);
  return {size, worth};
}


void DLVPartition::getGroupComp(vector<long long>& ids, int layer, long long group_id){
  string current_ptable = getPName(layer);
  char* vals[1];
  assign(vals, 0, group_id);
  _res = PQexecPrepared(_conn, fmt::format("comp_{}", current_ptable).c_str(), 1, vals, NULL, NULL, 0);
  free(vals, 1);
  for (int i = 0; i < PQntuples(_res); i ++) ids.push_back(atoll(PQgetvalue(_res, i, 0)));
  PQclear(_res);
}

vector<pair<double, double>> DLVPartition::getGroupIntervals(int layer, long long group_id){
  vector<pair<double, double>> intervals;
  string current_initial_gtable = getInitialGName(layer);
  char* vals[1];
  assign(vals, 0, group_id);
  _res = PQexecPrepared(_conn, fmt::format("interval_{}", current_initial_gtable).c_str(), 1, vals, NULL, NULL, 0);
  free(vals, 1);
  for (int i = 0; i < PQnfields(_res); i ++){
    intervals.push_back(atop(PQgetvalue(_res, 0, i)));
  }
  PQclear(_res);
  return intervals;
}

long long DLVPartition::getGroupContaining(int layer, vector<double> values){
  string current_initial_gtable = getInitialGName(layer);
  int m = (int) values.size();
  char* vals[m];
  for (int i = 0; i < m; i ++){
    assign(vals, i, values[i]);
  }
  _res = PQexecPrepared(_conn, fmt::format("group_{}", current_initial_gtable).c_str(), m, vals, NULL, NULL, 0);
  free(vals, m);
  long long group_id = -1;
  if (PQntuples(_res)) group_id = atoll(PQgetvalue(_res, 0, 0));
  PQclear(_res);
  return group_id;
}

void DLVPartition::getNeighboringGroupsRecurse(unordered_set<long long> &group_ids, vector<double> &values, int layer, int index, vector<pair<double, double>> &intervals){
  if (index == (int) intervals.size()){
    group_ids.insert(getGroupContaining(layer, values));
  } else{
    auto [left, right] = intervals[index];
    if (left == -DBL_MAX && right == DBL_MAX){
      values[index] = 0.0;
      getNeighboringGroupsRecurse(group_ids, values, layer, index+1, intervals);
    } else if (left == -DBL_MAX){
      values[index] = right;
      getNeighboringGroupsRecurse(group_ids, values, layer, index+1, intervals);
      values[index] = right + kIntervalEps;
      getNeighboringGroupsRecurse(group_ids, values, layer, index+1, intervals);
    } else if (right == DBL_MAX){
      values[index] = left;
      getNeighboringGroupsRecurse(group_ids, values, layer, index+1, intervals);
      values[index] = left - kIntervalEps;
      getNeighboringGroupsRecurse(group_ids, values, layer, index+1, intervals);
    } else{
      values[index] = (left + right) / 2;
      getNeighboringGroupsRecurse(group_ids, values, layer, index+1, intervals);
      values[index] = left - kIntervalEps;
      getNeighboringGroupsRecurse(group_ids, values, layer, index+1, intervals);
      values[index] = right + kIntervalEps;
      getNeighboringGroupsRecurse(group_ids, values, layer, index+1, intervals);
    }
  }
}

void DLVPartition::getNeighboringGroups(unordered_set<long long> &group_ids, int layer, long long group_id){
  vector<pair<double, double>> intervals = getGroupIntervals(layer, group_id);
  int neighbor_counts = 1;
  for (auto &p : intervals){
    int multiplier = 3;
    if (p.first == -DBL_MAX) multiplier --;
    if (p.second == DBL_MAX) multiplier --;
    neighbor_counts *= multiplier;
  }
  group_ids.reserve(neighbor_counts);
  vector<double> values (cols.size());
  getNeighboringGroupsRecurse(group_ids, values, layer, 0, intervals);
  group_ids.erase(group_id);
}