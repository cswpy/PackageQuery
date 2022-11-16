#include "lsr.h"

#include "pb/core/dual_reducer.h"

const int kSleepPeriod = 25; // In ms

void LayeredSketchRefine::init(){
  pg = new PgManager();
  _conn = PQconnectdb(pg->conninfo.c_str());
  assert(PQstatus(_conn) == CONNECTION_OK);
  _res = NULL;
  status = NotFound;
  INIT_CLOCK(pro);
}

DLVPartition* LayeredSketchRefine::getDLVPartition(LsrProb* prob){
  _sql = fmt::format("SELECT cols, group_ratio, lp_size, layer_count FROM {} WHERE table_name='{}' AND partition_name='{}';", kPartitionTable, prob->det_sql.table_name, prob->partition_name);
  _res = PQexec(_conn, _sql.c_str());
  if (PQntuples(_res)){
    DLVPartition* partition = new DLVPartition(prob, pgStringSplit(PQgetvalue(_res, 0, 0)), 
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

void LayeredSketchRefine::formulateDetProb(int core, LsrProb &prob, DetProb &det_prob, string current_gtable, const vector<long long> &ids){
  int m = prob.det_sql.att_cols.size() + 1;
  int n = (int) ids.size();
  det_prob.resize(m, n);
  det_prob.u.fill((double) prob.det_sql.u);
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
  string col_names = join(prob.det_sql.att_cols, ",");
  #pragma omp parallel num_threads(core)
  {
    CREATE_CLOCK(local_pro);
    string sql;
    PGconn* conn = PQconnectdb(pg->conninfo.c_str());
    assert(PQstatus(conn) == CONNECTION_OK);
    PGresult *res = NULL;
    int seg = omp_get_thread_num();
    int start_id = seg * chunk;
    int end_id = min((seg + 1) * chunk - 1, n - 1);
    string string_ids = "";
    for (int i = start_id; i < end_id; i ++) string_ids += to_string(ids[i]) + ",";
    if (end_id >= start_id) {
      string_ids += to_string(ids[end_id]);
      sql = fmt::format("SELECT {},{},{} FROM \"{}\" WHERE {} IN ({});",
        kId, prob.det_sql.obj_col, col_names, current_gtable, kId, string_ids);
      START_CLOCK(local_pro, 0);
      res = PQexec(conn, sql.c_str());
      END_CLOCK(local_pro, 0);
      for (int i = 0; i < PQntuples(res); i++){
        int index = i + start_id;
        det_prob.ids[index] = atoll(PQgetvalue(res, i, 0));
        det_prob.c(index) = atof(PQgetvalue(res, i, 1));
        for (int j = 0; j < m-1; j ++){
          det_prob.A(j, index) = atof(PQgetvalue(res, i, 2+j)); 
        }
        det_prob.A(m-1, index) = 1.0;
      }
      ADD_CLOCK(local_pro, core);
    }
    PQclear(res);
    PQfinish(conn);
  }
  det_prob.truncate();
}

LayeredSketchRefine::LayeredSketchRefine(int core, LsrProb &prob, bool is_safe){
  init();
  std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
  partition = getDLVPartition(&prob);
  if (!partition){
    status = NoPartitionFound;
    return;
  }
  if (!partition->isCompatible()){
    status = IncompatiblePartition;
    return;
  }

  double limit_size_per_group = DBL_MAX;
  if (isLess(kOutlierPercentage, 1.0)){
    limit_size_per_group = 
      (1 - kOutlierPercentage*partition->group_ratio)
      /(partition->group_ratio*(1-kOutlierPercentage));
  }

  // Phase-0: Filtering
  string current_gtable = partition->getGName(partition->layer_count);
  int n = pg->getSize(current_gtable);
  vector<long long> last_layer_ids;
  if (partition->is_filtering){
    // Filtering
    vector<string> att_names;
    for (int i = 0; i < (int) partition->query_cols.size(); i ++){
      att_names[i] = fmt::format("{} DOUBLE PRECISION", partition->query_cols[i]);
    }
    #if DEBUG
      MeanVar global_mv = MeanVar(partition->query_cols.size());
    #endif
    string gname;
    long long group_count, chunk;
    int start_index = 0;
    #pragma omp parallel num_threads(core)
    {
      string sql;
      PGconn *conn = PQconnectdb(pg->conninfo.c_str());
      assert(PQstatus(conn) == CONNECTION_OK);
      PGresult *res = NULL;
      DLVPartition *loc_partition = new DLVPartition(&prob, partition->cols, partition->group_ratio, partition->lp_size, partition->layer_count);
      
      int local_start_index = 0;
      vector<long long> local_ids;
      local_ids.reserve(n);

      for (int layer = 1; layer <= partition->layer_count; layer ++){
        #pragma omp master
        {
          gname = partition->getGName(layer);
          pg->dropTable(gname);
          _sql = fmt::format(""
            "CREATE TABLE IF NOT EXISTS \"{}\"("
            "	{} BIGINT,"
            " {}"
            ");", gname, kId, join(att_names, ","));
          _res = PQexec(_conn, _sql.c_str());
          PQclear(_res);

          group_count = pg->getSize(partition->getInitialGName(layer));
          chunk = ceilDiv(group_count, (long long) core);
        }
        #pragma omp barrier
        int seg = omp_get_thread_num();
        long long start_id = seg * chunk + 1;
        long long end_id = min((seg + 1) * chunk, group_count);
        if (end_id >= start_id){
          sql = fmt::format("COPY \"{}\" FROM STDIN with(delimiter ',');", gname);
          res = PQexec(conn, sql.c_str());
          assert(PQresultStatus(res) == PGRES_COPY_IN);
          PQclear(res);
          if (layer == 1){
            // Fast filtering if the intervals lie outside the filtering conditions
            unordered_map<string, int> cols_map;
            for (int i = 0; i < (int) partition->cols.size(); i ++){
              cols_map[partition->cols[i]] = i;
            }
            for (long long group_id = start_id; group_id <= end_id; group_id ++){
              vector<pair<double, double>> intervals = loc_partition->getGroupIntervals(layer, group_id);
              bool is_filtered = false;
              for (int j = 0; j < (int) prob.det_sql.filter_cols.size(); j++){
                int col_index = cols_map[prob.det_sql.filter_cols[j]];
                if (!doesIntersect(intervals[col_index], prob.det_sql.filter_intervals[j])){
                  is_filtered = true;
                  break;
                }
              }
              if (!is_filtered){
                // Go through group and check individual original tuples
                MeanVar mv = MeanVar(loc_partition->query_cols.size());
                loc_partition->getGroupFilteredStat(mv, layer, group_id);
                if (mv.sample_count){
                  string data = fmt::format("{},{}\n", group_id, join(mv.getVar(), kPrecision));
                  assert(PQputCopyData(_conn, data.c_str(), (int) data.length()) == 1);
                }
                #if DEBUG
                #pragma omp critical (c6)
                {
                  global_mv.add(mv);
                }
                #endif
                if (layer == partition->layer_count){
                  local_ids.push_back(group_id);
                }
              }
            }
          } else{
            // Normal filtering since the intervals are no longer correct for original tuples
            for (long long group_id = start_id; group_id <= end_id; group_id ++){
              MeanVar mv = MeanVar(loc_partition->query_cols.size());
              loc_partition->getGroupFilteredStat(mv, layer, group_id);
              if (mv.sample_count){
                string data = fmt::format("{},{}\n", group_id, join(mv.getVar(), kPrecision));
                assert(PQputCopyData(_conn, data.c_str(), (int) data.length()) == 1);
              }
              if (layer == partition->layer_count){
                local_ids.push_back(group_id);
              }
            }
          }
          assert(PQputCopyEnd(_conn, NULL) == 1);
          _res = PQgetResult(_conn);
          assert(PQresultStatus(_res) == PGRES_COMMAND_OK);
          PQclear(_res);
        }
        #pragma omp barrier
        #pragma omp master
        {
          _sql = fmt::format("CREATE INDEX \"{}_id\" ON \"{}\" USING btree ({});", gname, gname, kId);
          _res = PQexec(_conn, _sql.c_str());
          PQclear(_res);
        }
        #pragma omp barrier
        if (layer < partition->layer_count) loc_partition->prepareGroupFilteredStatSql(layer+1);
        #pragma omp barrier
      }
      PQfinish(conn);
      delete loc_partition;
      #pragma omp critical (c7)
      {
        local_start_index = start_index;
        start_index += local_ids.size();
      }
      #pragma omp barrier
      #pragma omp master
      {
        last_layer_ids.resize(start_index);
      }
      #pragma omp barrier
      if (local_ids.size()){
        memcpy(&(last_layer_ids[local_start_index]), &(local_ids[0]), local_ids.size()*sizeof(long long));
      }
    }
    #if DEBUG
      VectorXd means = global_mv.getMean();
      means.conservativeResize(prob.det_sql.att_cols.size());
      VectorXd vars = global_mv.getVar();
      vars.conservativeResize(prob.det_sql.att_cols.size());
      prob.rebound(means, vars);
    #endif
  } else{
    last_layer_ids.resize(n);
    iota(last_layer_ids.begin(), last_layer_ids.end(), 1);
  }
  DetProb det_prob;
  // Phase-1: Retrieve last layer
  formulateDetProb(core, prob, det_prob, current_gtable, last_layer_ids);
  // Phase-2: Iteratively go to next layer
  for (int layer = partition->layer_count; layer >= 1; layer --){
    int n = (int) det_prob.c.size();
    // Phase-2a: Sketch
    DualReducer dr = DualReducer(core, det_prob, is_safe);
    if (dr.status != Found){
      status = dr.status;
      return;
    }
    priority_queue<pair<double, long long>> pq;
    unordered_set<long long> group_ids_set;
    vector<long long> ids;
    long long total_size = 0;
    long long global_start_index = 0;
    int active_count = 0;
    #pragma omp parallel num_threads(core)
    {
      DLVPartition *loc_partition = new DLVPartition(&prob, partition->cols, partition->group_ratio, partition->lp_size, partition->layer_count);
      vector<long long> loc_ids; loc_ids.reserve(loc_partition->lp_size);
      #pragma omp for
      for (int i = 0; i < n; i ++){
        // Condition for sketch
        if (dr.lp_sol(i) > 0 || dr.ilp_sol(i) > 0){
          long long group_id = det_prob.ids[i];
          if (total_size <= partition->lp_size){
            long long size = loc_partition->getGroupComp(loc_ids, layer, group_id, limit_size_per_group);
            #pragma omp atomic
            total_size += size;
            #pragma omp critical (c2)
            {
              pq.emplace(det_prob.c(i), group_id);
              group_ids_set.insert(group_id);
            }
          }
        }
      }
      // Phase-2b: Shade
      bool is_active = false;
      while (true){
        if (total_size > partition->lp_size) break;
        long long group_id = -1;
        #pragma omp critical (c3)
        {
          if (pq.size()){
            group_id = pq.top().second;
            pq.pop();
            if (!is_active){
              active_count ++;
              // cout << "CURRENT ACTIVE: " << active_count << endl;
              is_active = true;
            }
          } else if (is_active){
            active_count --;
            is_active = false;
          }
        }
        if (group_id < 0){
          // is_active = false
          std::this_thread::sleep_for(std::chrono::milliseconds(kSleepPeriod));
          if (!active_count) break;
        } else{
          unordered_set<long long> group_ids;
          loc_partition->getNeighboringGroups(group_ids, layer, group_id);
          for (const auto& gid : group_ids){
            if (total_size <= partition->lp_size){
              bool is_adding = false;
              #pragma omp critical (c4)
              {
                if (group_ids_set.find(gid) == group_ids_set.end()){
                  is_adding = true;
                  group_ids_set.insert(gid);
                }
              }
              if (is_adding){
                double worth = loc_partition->getGroupWorthness(layer, gid);
                long long actual_size = loc_partition->getGroupComp(loc_ids, layer, gid, limit_size_per_group);
                #pragma omp atomic
                total_size += actual_size;
                #pragma omp critical (c3)
                {
                  pq.emplace(worth, gid);
                }
              }
            }
          }
          // Debug
          // #pragma omp critical
          // {
          //   cout << "ORIGIN: " << group_id << endl;
          //   for (const auto &elem : group_ids) cout << elem << " ";
          //   cout << endl;
          // }
        }
      }

      // Debug
      // #pragma omp critical
      // {
      //   cout << "SIZE: " << loc_ids.size() << endl;
      //   for (auto id : loc_ids){
      //     cout << id << " ";
      //   }
      //   cout << endl;
      // }

      // Phase-2c: Refine
      long long start_index = 0;
      #pragma omp critical (c5)
      {
        start_index = global_start_index;
        global_start_index += loc_ids.size();
      }
      #pragma omp barrier
      #pragma omp master
      {
        ids.resize(total_size);
      }
      #pragma omp barrier
      if (loc_ids.size()) memcpy(&(ids[start_index]), &(loc_ids[0]), loc_ids.size()*sizeof(long long));
      ADD_CLOCK(loc_partition->pro, core);
      delete loc_partition;
    }

    string next_gtable;
    if (layer > 1) next_gtable = partition->getGName(layer-1);
    else next_gtable = prob.det_sql.table_name;
    formulateDetProb(core, prob, det_prob, next_gtable, ids);
    // Debug
    // for (auto v : ids) cout << v << " ";
    // cout << endl;
    // cout << ids.size() << endl;
  }
  if (partition->is_filtering){
    for (int layer = 1; layer <= partition->layer_count; layer ++){
      string gname = partition->getGName(layer);
      pg->dropTable(gname);
    }
  }
  // Phase-3: Final solution
  DualReducer dr = DualReducer(core, det_prob, is_safe);
  if (dr.status != Found){
    status = dr.status;
    return;
  }
  status = Found;
  ilp_score = lp_score = 0;
  for (int i = 0; i < (int) det_prob.c.size(); i ++){
    if (isGreater(dr.ilp_sol(i), 0)){
      ilp_sol[det_prob.ids[i]] = (long long) dr.ilp_sol(i);
      ilp_score += ilp_sol[det_prob.ids[i]]*det_prob.c(i);
    }
    if (isGreater(dr.lp_sol(i), 0)){
      lp_sol[det_prob.ids[i]] = dr.lp_sol(i);
      lp_score += lp_sol[det_prob.ids[i]]*det_prob.c(i);
    }
  }
  exe_ilp = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start).count() / 1000000.0;
  exe_lp = exe_ilp - (dr.exe_ilp - dr.exe_lp);
}