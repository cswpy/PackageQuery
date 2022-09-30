#include "lsr.h"

#include "pb/util/uconfig.h"
#include "pb/core/dual_reducer.h"

const string kTempPrefix = "tmp";

DLVPartition::~DLVPartition(){
  delete pg;
  PQfinish(_conn);
}

DLVPartition::DLVPartition(string table_name, string partition_name, vector<string> cols, double group_ratio, long long lp_size, int layer_count)
  : table_name(table_name), partition_name(partition_name), cols(cols), group_ratio(group_ratio), lp_size(lp_size), layer_count(layer_count){
  pg = new PgManager();
  _conn = PQconnectdb(pg->conninfo.c_str());
  assert(PQstatus(_conn) == CONNECTION_OK);
  _res = NULL;
}

string DLVPartition::getPName(int layer){
  return fmt::format("[{}P]{}_{}", layer, table_name, partition_name);
}

string DLVPartition::getGName(int layer, bool is_filtered){
  if (!is_filtered) return fmt::format("[{}G]{}_{}", layer, table_name, partition_name);
  else return fmt::format("{}_[{}G]{}_{}", kTempPrefix, layer, table_name, partition_name);
}

bool DLVPartition::isCompatible(const LsrProb &prob){
  bool res = isIn(cols, prob.obj_col);
  for (string s : prob.cols) res &= isIn(cols, s);
  vector<string> col_names = pg->listColumns(prob.table_name);
  for (string s : prob.filter_cols) res &= isIn(col_names, s);
  return res;
}

long long DLVPartition::getGroupSize(int layer, long long id, bool is_filtered){
  string current_gtable = getGName(layer, is_filtered);
  return 0;
}

void LayeredSketchRefine::init(){
  vector<string> names = {"Init", "All", "All", "FetchData", "ProcessData", "WriteData", "CreateIndex", "CreateTable", "gist", "id", "tid", "gid"};
  pro = Profiler(names);
  pro.clock(0);
  pg = new PgManager();
  _conn = PQconnectdb(pg->conninfo.c_str());
  assert(PQstatus(_conn) == CONNECTION_OK);
  _res = NULL;
  status = NotFound;
  pro.stop(0);
}

DLVPartition* LayeredSketchRefine::getDLVPartition(string table_name, string partition_name){
  _sql = fmt::format("SELECT cols, group_ratio, lp_size, layer_count FROM {} WHERE table_name='{}' AND partition_name='{}';", kPartitionTable, table_name, partition_name);
  _res = PQexec(_conn, _sql.c_str());
  if (PQntuples(_res)){
    DLVPartition* partition = new DLVPartition(table_name, partition_name, pgStringSplit(PQgetvalue(_res, 0, 0)), 
      atof(PQgetvalue(_res, 0, 1)), atoll(PQgetvalue(_res, 0, 2)), atoi(PQgetvalue(_res, 0, 3)));
    PQclear(_res);
    return partition;
  }
  return NULL;
}

LayeredSketchRefine::~LayeredSketchRefine(){
  PQfinish(_conn);
  if (partition) delete partition;
  delete pg;
}

void LayeredSketchRefine::formulateDetProb(int core, const LsrProb &prob, DetProb &det_prob, string current_gtable, const vector<long long> &ids){
  int m = prob.cols.size() + 1;
  int n = (int) ids.size();
  det_prob.resize(m, n);
  det_prob.u.fill((double) prob.u);
  for (int i = 0; i < m; i ++){
    if (i < m - 1){
      det_prob.bl(i) = prob.bl(i);
      det_prob.bu(i) = prob.bu(i);
    } else{
      det_prob.bl(i) = prob.cl;
      det_prob.bu(i) = prob.cu;
    }
  }
  int chunk = ceilDiv(n, core);
  string col_names = join(prob.cols, ",");
  #pragma omp parallel num_threads(core)
  {
    string sql;
    PGconn* conn = PQconnectdb(pg->conninfo.c_str());
    assert(PQstatus(conn) == CONNECTION_OK);
    PGresult *res = NULL;
    int seg = omp_get_thread_num();
    int start_id = seg * chunk;
    int end_id = min((seg + 1) * chunk - 1, n - 1);
    string string_ids = "";
    for (int i = start_id; i < end_id; i ++) string_ids += to_string(ids[i]) + ",";
    if (end_id >= start_id) string_ids += to_string(ids[end_id]);
    sql = fmt::format("SELECT {},{},{} FROM \"{}\" WHERE {} IN ({});",
      kId, prob.obj_col, col_names, current_gtable, kId, string_ids);
    res = PQexec(conn, sql.c_str());
    for (int i = 0; i < PQntuples(res); i++){
      int index = i + start_id;
      det_prob.ids[index] = atoll(PQgetvalue(res, i, 0));
      det_prob.c(index) = atof(PQgetvalue(res, i, 1));
      for (int j = 0; j < m-1; j ++){
        det_prob.A(j, index) = atof(PQgetvalue(res, i, 2+j)); 
      }
      det_prob.A(m-1, index) = 1.0;
    }
    PQclear(res);
    PQfinish(conn);
  }
  det_prob.truncate();
}

LayeredSketchRefine::LayeredSketchRefine(int core, const LsrProb &prob){
  init();
  pro.clock(1);
  partition = getDLVPartition(prob.table_name, prob.partition_name);
  if (!partition){
    status = NoPartitionFound;
    return;
  }
  if (!partition->isCompatible(prob)){
    status = IncompatiblePartition;
    return;
  }
  // Phase-0: Filtering
  bool is_filtering = prob.filter_cols.size() > 0;
  if (is_filtering){
    // Filtering
  }
  DetProb det_prob;
  if (is_filtering){
  } else{ // No filter
    // Phase-1: Retrieve last layer
    {
      string current_gtable = partition->getGName(partition->layer_count);
      int n = pg->getSize(current_gtable);
      vector<long long> ids (n);
      iota(ids.begin(), ids.end(), 1);
      formulateDetProb(core, prob, det_prob, current_gtable, ids);
    }
    // Phase-2: Iteratively go to next layer
    for (int layer = partition->layer_count; layer >= 1; layer --){
      int m = (int) det_prob.bl.size();
      int n = (int) det_prob.c.size();
      // Phase-2a: Sketch
      DualReducer dr = DualReducer(core, det_prob);
      if (dr.status != Found){
        status = dr.status;
        return;
      }
      priority_queue<pair<double, long long>> pq;
      unordered_set<long long> ids_set;
      long long total_size = 0;
      #pragma omp parallel num_threads(core)
      {
        DLVPartition *loc_partition = new DLVPartition(partition->table_name, partition->partition_name, partition->cols, partition->group_ratio, partition->lp_size, partition->layer_count);
        #pragma omp for
        for (int i = 0; i < n; i ++){
          // Condition for sketch
          if (dr.lp_sol(i) > 0 || dr.ilp_sol(i) > 0){
            #pragma omp critical
            {
              pq.emplace(det_prob.c(i), det_prob.ids[i]);
              ids_set.insert(det_prob.ids[i]);
            }
          }
        }
        delete loc_partition;
      }
      // Phase-2b: Shade
      break;
    }
  }
  pro.stop(1);
}