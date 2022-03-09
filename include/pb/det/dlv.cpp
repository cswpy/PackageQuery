#define FMT_HEADER_ONLY

#include "fmt/core.h"
#include "dlv.h"
#include "pb/util/uconfig.h"
#include "pb/util/umisc.h"
#include "pb/util/unumeric.h"
#include "pb/util/udeclare.h"
#include "pb/core/map_sort.h"

using namespace pb;

const string kIntervalType = "floatrange";
const string kTempTable = "tmp";
const double eps = 1e-16;
const int kTempReserveSize = 1000;
// const string kTraverseFunction = "traverse";
// const string kClearLockFunction = "clear_lock";
// const string kProcessReserveFunction = "process_reserve";
// const int kReserveSize = (int) ceilDiv(kInMemorySize, kPCore);

#define at(row, col) (m*(col)+(row))

struct IndexComp{
  const double* mat;
  int j, m;
  IndexComp(const double* mat, int j, int m): mat(mat), j(j), m(m){
  }
  inline bool operator()(int i1, int i2){
    return mat[at(j, i1)] < mat[at(j, i2)];
  }
};

double DynamicLowVariance::kGroupRatio = 0.01;
double DynamicLowVariance::kVarScale = 2.5;
int DynamicLowVariance::kMaxSize = 100000;
// int DynamicLowVariance::kInitialSize = (int) ceil(kPCore / (kGroupRatio * kGroupRatio));

DynamicLowVariance::~DynamicLowVariance(){
  PQfinish(_conn);
  delete pg;
}

DynamicLowVariance::DynamicLowVariance(string dbname): dbname(dbname){
  vector<string> names = {"Init", "ComputeStat", "PartialPartition", "FetchData", "ProcessData", "WriteData", "CreateIndex", "CreateTable", "gist", "id", "tid", "gid"};
  pro = Profiler(names);
  init();
}

void DynamicLowVariance::init(){
  pro.clock(0);
  pg = new PgManager(dbname);

  string conninfo = fmt::format("postgresql://{}@{}?port={}&dbname={}&password={}", kPgUser, kPgHostaddr, kPgPort, dbname, kPgPassword);
  _conn = PQconnectdb(conninfo.c_str());
  assert(PQstatus(_conn) == CONNECTION_OK);
  _res = NULL;

  _sql = fmt::format("SET client_min_messages = warning;");
  _res = PQexec(_conn, _sql.c_str());
  PQclear(_res);

  _sql = fmt::format(""\
    "DO $$ BEGIN "\
      "CREATE TYPE {} AS RANGE ("\
      "subtype = float8,"\
      "subtype_diff = float8mi"\
      ");"
    "EXCEPTION "\
      " WHEN duplicate_object THEN null;"\
    "END $$;", kIntervalType);
  _res = PQexec(_conn, _sql.c_str());
  PQclear(_res);

  // _sql = fmt::format(""\
  //   "CREATE OR REPLACE FUNCTION {}() "\
  //   "RETURNS void "\
  //   "LANGUAGE plpgsql AS "\
  //   "$$ "\
  //   "DECLARE "\
  //   "	rec record; "\
  //   "BEGIN "\
  //   "	FOR rec IN SELECT pid FROM pg_locks WHERE locktype='advisory' LOOP "\
  //   "		PERFORM pg_terminate_backend(rec.pid); "\
  //   "	END LOOP; "\
  //   "END; "\
  //   "$$; ", kClearLockFunction);
  // _res = PQexec(_conn, _sql.c_str());
  // PQclear(_res);

  _sql = fmt::format(""
    "CREATE TABLE IF NOT EXISTS {}("
    "	table_name VARCHAR(63) UNIQUE NOT NULL,"
    "	size BIGINT,"
    " cols TEXT[],"
    "	mean DOUBLE PRECISION[],"
    "	M2 DOUBLE PRECISION[],"
    "	amin DOUBLE PRECISION[],"
    "	amax DOUBLE PRECISION[]"
    ")", kStatTable);
  _res = PQexec(_conn, _sql.c_str());
  PQclear(_res);

  _sql = fmt::format(""
    "CREATE TABLE IF NOT EXISTS {}("
    "	table_name VARCHAR(63) NOT NULL,"
    "	partition_name VARCHAR(31) NOT NULL,"
    "	cols TEXT[],"
    "	group_ratio DOUBLE PRECISION,"
    "	max_size INTEGER,"
    " layer_count INTEGER,"
    " UNIQUE (table_name, partition_name)"
    ");", kPartitionTable);
  _res = PQexec(_conn, _sql.c_str());
  PQclear(_res);
  pro.stop(0);
}

void DynamicLowVariance::writeStats(string table_name, Stat *stat){
  _sql = fmt::format(""
    "INSERT INTO {}(table_name, size, cols, mean, M2, amin, amax) "
    "VALUES ('{}', {}, {}, {}, {}, {}, {}) "
    "ON CONFLICT (table_name) "
    "DO UPDATE SET size=EXCLUDED.size,cols=EXCLUDED.cols,mean=EXCLUDED.mean,M2=EXCLUDED.M2,amin=EXCLUDED.amin,amax=EXCLUDED.amax;", 
    kStatTable, table_name, stat->size, pgJoin(stat->cols), pgJoin(stat->mean, kPrecision), pgJoin(stat->M2, kPrecision), pgJoin(stat->amin, kPrecision), pgJoin(stat->amax, kPrecision));
  _res = PQexec(_conn, _sql.c_str());
  PQclear(_res);
}

Stat* DynamicLowVariance::readStats(string table_name){
  _sql = fmt::format("SELECT size, cols, mean, M2, amin, amax FROM \"{}\" WHERE table_name='{}';", kStatTable, table_name);
  _res = PQexec(_conn, _sql.c_str());
  Stat* stat = new Stat(pgStringSplit(PQgetvalue(_res, 0, 1)));
  stat->add(atoll(PQgetvalue(_res, 0, 0)), pgValueSplit(PQgetvalue(_res, 0, 2)), pgValueSplit(PQgetvalue(_res, 0, 3)), pgValueSplit(PQgetvalue(_res, 0, 4)), pgValueSplit(PQgetvalue(_res, 0, 5)));
  PQclear(_res);
  return stat;
}

void DynamicLowVariance::partition(string table_name, string partition_name){
  pro.clock(1);
  vector<string> cols = pg->getNumericCols(table_name);
  // Stat *stat = pg->computeStats(table_name, cols);
  // writeStats(table_name, stat);
  // delete stat;
  pro.stop(1);
  pro.clock(2);
  partition(table_name, partition_name, cols);
  pro.stop(2);
}

void DynamicLowVariance::partition(string table_name, string partition_name, const vector<string> &cols){
  doPartition(table_name, partition_name, cols);
}

void DynamicLowVariance::quickPartition(string table_name, Stat *stat, const vector<string> &cols){
  // Assumes uniform distribtion
  int max_var_index = -1;
  double max_var = -1;
  for (int i = 0; i < cols.size(); i ++){
    double var = stat->getVar(cols[i]);
    if (max_var < var){
      max_var = var;
      max_var_index = i;
    }
  }
  double min_att = stat->amin(max_var_index);
  double max_att = stat->amax(max_var_index);

  long long size = stat->size;
  long long chunk = ceilDiv(size, kPCore);
  long long bucket = ceilDiv(size, kInMemorySize);
  {
    string sql;
    string conninfo = fmt::format("postgresql://{}@{}?port={}&dbname={}&password={}", kPgUser, kPgHostaddr, kPgPort, dbname, kPgPassword);
    PGconn* conn = PQconnectdb(conninfo.c_str());
    assert(PQstatus(conn) == CONNECTION_OK);
    PGresult *res = NULL;

    for (long long i = 0; i < bucket; i ++){
      sql = fmt::format(""
        "CREATE TABLE IF NOT EXISTS {}{}("
        "	tid BIGINT"
        ");", kTempTable, i);
      res = PQexec(conn, sql.c_str());
      PQclear(res);
    }
    PQfinish(conn);
  }
  #pragma omp parallel num_threads(kPCore)
  {
    string sql;
    string conninfo = fmt::format("postgresql://{}@{}?port={}&dbname={}&password={}", kPgUser, kPgHostaddr, kPgPort, dbname, kPgPassword);
    PGconn* conn = PQconnectdb(conninfo.c_str());
    assert(PQstatus(conn) == CONNECTION_OK);
    PGresult *res = NULL;

    PGconn* _conn = PQconnectdb(conninfo.c_str());
    assert(PQstatus(_conn) == CONNECTION_OK);
    PGresult *_res = NULL;

    int seg = omp_get_thread_num();
    long long start_id = seg * chunk + 1;
    long long end_id = min((seg + 1) * chunk, size);
    vector<int> sz (bucket, 0);
    vector<vector<long long>> cache (bucket, vector<long long>(kTempReserveSize));

    sql = fmt::format("SELECT {},{} FROM \"{}\" WHERE {} BETWEEN {} AND {};", kId, stat->cols[max_var_index], table_name, kId, start_id, end_id);
    res = PQexec(conn, sql.c_str());
    for (int i = 0; i < PQntuples(res); i++){
      long long tid = atol(PQgetvalue(res, i, 0));
      long long bucket_ind = min((long long) floor((atof(PQgetvalue(res, i, 1)) - min_att) / (max_att - min_att + eps) * bucket), bucket-1);
      cache[bucket_ind][sz[bucket_ind]] = tid;
      sz[bucket_ind] ++;
      if (sz[bucket_ind] == kTempReserveSize){
        sql = fmt::format("COPY {}{} FROM STDIN with(delimiter ',');", kTempTable, bucket_ind);
        _res = PQexec(_conn, sql.c_str());
        assert(PQresultStatus(_res) == PGRES_COPY_IN);
        PQclear(_res);
        string data = "";
        for (int j = 0; j < sz[bucket_ind]; j ++) data += fmt::format("{}\n", cache[bucket_ind][j]);
        assert(PQputCopyData(_conn, data.c_str(), (int) data.length()) == 1);
        assert(PQputCopyEnd(_conn, NULL) == 1);
        _res = PQgetResult(_conn);
        assert(PQresultStatus(_res) == PGRES_COMMAND_OK);
        sz[bucket_ind] = 0;
      }
    }
    PQclear(res);

    for (long long i = 0; i < bucket; i ++){
      if (sz[i]){
        sql = fmt::format("COPY {}{} FROM STDIN with(delimiter ',');", kTempTable, i);
        _res = PQexec(_conn, sql.c_str());
        assert(PQresultStatus(_res) == PGRES_COPY_IN);
        PQclear(_res);
        string data = "";
        for (int j = 0; j < sz[i]; j ++) data += fmt::format("{}\n", cache[i][j]);
        assert(PQputCopyData(_conn, data.c_str(), (int) data.length()) == 1);
        assert(PQputCopyEnd(_conn, NULL) == 1);
        _res = PQgetResult(_conn);
        assert(PQresultStatus(_res) == PGRES_COMMAND_OK);
      }
    }
    
    PQfinish(conn);
    PQfinish(_conn);
  }

  {
    string sql;
    string conninfo = fmt::format("postgresql://{}@{}?port={}&dbname={}&password={}", kPgUser, kPgHostaddr, kPgPort, dbname, kPgPassword);
    PGconn* conn = PQconnectdb(conninfo.c_str());
    assert(PQstatus(conn) == CONNECTION_OK);
    PGresult *res = NULL;

    for (long long i = 0; i < bucket; i ++){
      sql = fmt::format("SELECT COUNT(*) FROM {}{};", kTempTable, i);
      res = PQexec(conn, sql.c_str());
      long long sz = atol(PQgetvalue(res, 0, 0));
      PQclear(res);
      cout << sz;
      if (i == bucket-1) cout << "\n";
      else cout << " ";

      sql = fmt::format("DROP TABLE IF EXISTS {}{};", kTempTable, i);
      res = PQexec(conn, sql.c_str());
      PQclear(res);
    }
    PQfinish(conn);
  }
}

void DynamicLowVariance::doPartition(string table_name, string suffix, const vector<string> &cols){
  Stat *stat = readStats(table_name);
  quickPartition(table_name, stat, cols);
  delete stat;
}

// void DynamicLowVariance::doPartition(string table_name, string suffix, const vector<string> &cols){
//   // Compute min_partition_var
//   int m = (int) cols.size();
//   Stat *stat = readStats(table_name);
//   long long n = stat->size;
//   if (n < kMaxSize) return;
//   long long total_groups = (long long) ceil(n * kGroupRatio);
//   double min_partition_var = 1.0;
//   for (string col : cols) min_partition_var *= stat->getVar(col);
//   delete stat;
//   min_partition_var /= total_groups;
//   min_partition_var /= total_groups;
//   min_partition_var = pow(min_partition_var, 1.0/m) * kVarScale * kVarScale;
//   int min_partition_size = (int) ceil(1.0 / (kGroupRatio * kGroupRatio));

//   int kInitialSize = n;
//   cout << "With " << kInitialSize << endl;
//   // Phase 1a
//   string symbolic_name = table_name;
//   if (suffix.length()) symbolic_name += "_" + suffix;
//   string p_name = nextPName(symbolic_name);
//   string g_name = nextGName(symbolic_name);
//   string col_names = join(cols, ",");
//   Stat st = Stat(cols);
//   double *A = new double [m*kInitialSize];
//   pair<double, int>* pis = new pair<double, int>[kInitialSize];
//   VectorXi group_indices (kInitialSize);
//   int partition_index = -1;
//   int chunk = (int) ceilDiv(kInitialSize, kPCore);
//   long long overall_chunk = ceilDiv(n, kPCore);
//   long long phase_two_chunk = ceilDiv(n - kInitialSize, kPCore);
//   pro.clock(3);
//   #pragma omp parallel num_threads(kPCore)
//   {
//     string conninfo = fmt::format("postgresql://{}@{}?port={}&dbname={}&password={}", kPgUser, kPgHostaddr, kPgPort, dbname, kPgPassword);
//     PGconn *conn = PQconnectdb(conninfo.c_str());
//     assert(PQstatus(conn) == CONNECTION_OK);
//     PGresult *res = NULL;
//     MeanVar mv = MeanVar(m);
//     int tn = omp_get_thread_num();
//     int left = tn * chunk;
//     int right = min((tn+1) * chunk, kInitialSize);
//     string select_sql = fmt::format("SELECT {} FROM \"{}\" WHERE {} BETWEEN {} AND {}", col_names, table_name, kId, left+1, right);
//     res = PQexec(conn, select_sql.c_str());
//     for (int i = 0; i < PQntuples(res); i ++){
//       int index = i + left;
//       for (int j = 0; j < m; j ++) A[at(j, index)] = atof(PQgetvalue(res, i, j));
//       mv.add(A+at(0, index));
//     }
//     PQclear(res);
//     PQfinish(conn);
//     #pragma omp critical
//     {
//       st.add(mv.sample_count, mv.getMean(), mv.getM2());
//     }
//     #pragma omp barrier
//     #pragma omp master
//     {
//       double max_total_var = -1;
//       for (int i = 0; i < m; i ++){
//         double var = st.getVar(i);
//         if (max_total_var < var){
//           max_total_var = var;
//           partition_index = i;
//         }
//       }
//     }
//     #pragma omp barrier
//     #pragma omp for nowait
//     for (int i = 0; i < kInitialSize; i ++) pis[i] = {A[at(partition_index, i)], i};
//   }
//   pro.stop(3);
//   pro.clock(4);
//   map_sort::Sort(pis, kInitialSize, kPCore);
//   int soft_group_lim = (int) ceil(kInitialSize * kGroupRatio);
//   double var_ratio = kGroupRatio * kGroupRatio * kVarScale;
//   int soft_partition_lim = (int) ceil(kVarScale / kGroupRatio);
//   int hard_group_lim = soft_group_lim + kPCore + soft_partition_lim;
//   long long group_count = kPCore;
//   int heap_length = kPCore;
//   vector<vector<double>> lows (m, vector<double>(hard_group_lim, -DBL_MAX));
//   vector<vector<double>> highs (m, vector<double>(hard_group_lim, DBL_MAX));
//   vector<tuple<double, int, int>> max_heap (hard_group_lim);
//   vector<VectorXi*> groups (hard_group_lim, NULL);

//   string condition_names, rec_names, t_names;
//   {
//     vector<string> conditions;
//     vector<string> recs;
//     vector<string> ts;
//     for (int i = 0; i < m; i ++){
//       conditions.push_back(fmt::format("t.interval_{} @> %L::float8", cols[i]));
//       recs.push_back(fmt::format("rec.{}", cols[i]));
//       ts.push_back(fmt::format("t.{}", cols[i]));
//     }
//     condition_names = join(conditions, " AND ");
//     rec_names = join(recs, ", ");
//     t_names = join(ts, ",");
//   }

//   #pragma omp parallel num_threads(kPCore)
//   {
//     {
//       int i = omp_get_thread_num();
//       int left = i * chunk;
//       int right = min((i+1)*chunk, kInitialSize);
//       groups[i] = new VectorXi(right - left);
//       if (i > 0) lows[partition_index][i] = A[at(partition_index, pis[left].second)];
//       if (i < kPCore - 1) highs[partition_index][i] = A[at(partition_index, pis[right].second)];

//       MeanVar mv = MeanVar(m);
//       for (int j = left; j < right; j ++){
//         int col_index = pis[j].second;
//         (*groups[i])(j-left) = col_index;
//         mv.add(A+at(0, col_index));
//       }
//       int local_pindex = -1;
//       VectorXd total_vars = mv.getM2();
//       double max_total_var = -1;
//       for (int j = 0; j < m; j ++){
//         double total_var = total_vars(j);
//         if (max_total_var < total_var){
//           max_total_var = total_var;
//           local_pindex = j;
//         }
//       }
//       max_heap[i] = {max_total_var, local_pindex, i};
//     }
//     #pragma omp barrier
//     #pragma omp master
//     {
//       delete[] pis;
//       make_heap(max_heap.begin(), max_heap.begin() + heap_length);
//     }
//     #pragma omp barrier
//     // Phase 1b
//     while (group_count < soft_group_lim){
//       int mi = -1;
//       int gi = -1; 
//       double max_total_var = 0;
//       #pragma omp critical
//       {
//         while (heap_length > 0 && max_total_var == 0){
//           tie(max_total_var, mi, gi) = max_heap[0];
//           pop_heap(max_heap.begin(), max_heap.begin() + heap_length);
//           heap_length --;
//         }
//       }
//       if (gi == -1) break;
//       auto& g = *groups[gi];
//       sort(g.begin(), g.end(), IndexComp(A, mi, m));
//       int delim_sz = soft_partition_lim;
//       vector<int> delims (delim_sz);
//       int delim_count = 0;
//       ScalarMeanVar smv = ScalarMeanVar();
//       double reduced_var = max_total_var / g.size() * var_ratio;
//       for (int i = 0; i < g.size(); i ++){
//         smv.add(A[at(mi, g(i))]);
//         if (smv.getVar() > reduced_var){
//           if (delim_count < delim_sz) delims[delim_count] = i;
//           else{
//             delim_sz += soft_partition_lim;
//             delims.resize(delim_sz);
//             delims[delim_count] = i;
//           }
//           delim_count ++;
//           smv.reset();
//           smv.add(A[at(mi, g(i))]);
//         }
//       }
//       if (delim_count < delim_sz) delims[delim_count] = (int) g.size();
//       else delims.push_back((int) g.size());
//       delim_count ++;

//       int g_start_index;
//       #pragma omp critical
//       {
//         g_start_index = (int) group_count;
//         group_count += delim_count-1;
//         if (group_count > hard_group_lim){
//           hard_group_lim += kPCore * soft_partition_lim;
//           for (int i = 0; i < m; i ++){
//             lows[i].resize(hard_group_lim, -DBL_MAX);
//             highs[i].resize(hard_group_lim, DBL_MAX);
//           }
//           max_heap.resize(hard_group_lim);
//           groups.resize(hard_group_lim, NULL);
//         }
//       }

//       for (int i = 0; i < delim_count-1; i ++){
//         int g_index = g_start_index + i;
//         int group_sz = delims[i+1] - delims[i];
//         VectorXi* gptr = new VectorXi(group_sz);
//         memcpy(&(*gptr)(0), &g(delims[i]), group_sz*sizeof(int));

//         groups[g_index] = gptr;
//         for (int j = 0; j < m; j ++){
//           lows[j][g_index] = lows[j][gi];
//           highs[j][g_index] = highs[j][gi];
//         }
//         lows[mi][g_index] = A[at(mi, g(delims[i]))];
//         if (i < delim_count-2) highs[mi][g_index] = A[at(mi, g(delims[i+1]))];
//         else highs[mi][g_index] = highs[mi][gi];

//         // Begin Repeatable code 
//         MeanVar mv = MeanVar(m);
//         for (int j = 0; j < group_sz; j ++){
//           mv.add(A+at(0, (*gptr)(j)));
//         }
//         int m_index = -1;
//         VectorXd total_vars = mv.getM2();
//         double m_var = 0;
//         for (int j = 0; j < m; j ++){
//           double total_var = total_vars(j);
//           if (m_var < total_var){
//             m_var = total_var;
//             m_index = j;
//           }
//         }
//         if (m_index != -1){
//           #pragma omp critical
//           {
//             max_heap[heap_length] = {m_var, m_index, g_index};
//             heap_length ++;
//             push_heap(max_heap.begin(), max_heap.begin() + heap_length);
//           }
//         }
//         // End Repeatable code 
//       }

//       {
//         // lows is good no need to change
//         if (delim_count != 1) highs[mi][gi] = A[at(mi, g(delims[0]))];
//         g.conservativeResize(delims[0]);

//         int group_sz = delims[0];
//         // Begin Repeatable code 
//         MeanVar mv = MeanVar(m);
//         for (int j = 0; j < group_sz; j ++){
//           mv.add(A+at(0, g(j)));
//         }
//         VectorXd total_vars = mv.getM2();
//         int m_index = -1;
//         double m_var = 0;
//         for (int j = 0; j < m; j ++){
//           double total_var = total_vars(j);
//           if (m_var < total_var){
//             m_var = total_var;
//             m_index = j;
//           }
//         }
//         if (m_index != -1){
//           #pragma omp critical
//           {
//             max_heap[heap_length] = {m_var, m_index, gi};
//             heap_length ++;
//             push_heap(max_heap.begin(), max_heap.begin() + heap_length);
//           }
//         }
//         // End Repeatable code 
//       }
//     }

//     // Phase 1c
//     #pragma omp master
//     {
//       pro.stop(4);
//       pro.clock(5);
//       pro.clock(7);
//       string conninfo = fmt::format("postgresql://{}@{}?port={}&dbname={}&password={}", kPgUser, kPgHostaddr, kPgPort, dbname, kPgPassword);
//       PGconn *conn = PQconnectdb(conninfo.c_str());
//       assert(PQstatus(conn) == CONNECTION_OK);
//       PGresult *res = NULL;

//       string drop_sql = fmt::format("DROP TABLE IF EXISTS \"{}\";", g_name);
//       res = PQexec(conn, drop_sql.c_str());
//       PQclear(res);

//       string atts_names = fmt::format("{} BIGINT,", kId);
//       for (auto col : cols){
//         atts_names += fmt::format(""\
//         "interval_{} {},"\
//         "{} DOUBLE PRECISION,", col, kIntervalType, col);
//       }
//       atts_names += "max_total_var DOUBLE PRECISION, max_size BIGINT, stsatus INTEGER";
//       string create_sql = fmt::format("CREATE TABLE IF NOT EXISTS \"{}\"({});", g_name, atts_names);
//       res = PQexec(conn, create_sql.c_str());
//       PQclear(res);

//       drop_sql = fmt::format("DROP TABLE IF EXISTS \"{}\";", p_name);
//       res = PQexec(conn, drop_sql.c_str());
//       PQclear(res);

//       create_sql = fmt::format("CREATE TABLE IF NOT EXISTS \"{}\"("\
//         "tid BIGINT,"\
//         "gid BIGINT"\
//         ");", p_name);
//       res = PQexec(conn, create_sql.c_str());
//       PQclear(res);

//       PQfinish(conn);
//       pro.stop(7);
//     }
//     #pragma omp barrier

//     {
//       string conninfo = fmt::format("postgresql://{}@{}?port={}&dbname={}&password={}", kPgUser, kPgHostaddr, kPgPort, dbname, kPgPassword);
//       PGconn *conn = PQconnectdb(conninfo.c_str());
//       assert(PQstatus(conn) == CONNECTION_OK);
//       PGresult *res = NULL;

//       // Populate G table
//       string sql = fmt::format("COPY \"{}\" FROM STDIN with(delimiter '|', null '{}');", g_name, kNullLiteral);
//       res = PQexec(conn, sql.c_str());
//       assert(PQresultStatus(res) == PGRES_COPY_IN);
//       PQclear(res);
//       #pragma omp for nowait
//       for (int i = 0; i < group_count; i ++){
//         const auto& g = *groups[i];
//         MeanVar mv = MeanVar(m);
//         string data = fmt::format("{}|", i+1);
//         for (int j : g) mv.add(A+at(0, j));
//         VectorXd mean = mv.getMean();
//         for (int j = 0; j < m; j ++){
//           data += fmt::format("[{:.{}Lf},{:.{}Lf}]|{:.{}Lf}|", lows[j][i], kPrecision, highs[j][i], kPrecision, mean(j), kPrecision);
//         }
//         data += fmt::format("{:.{}Lf}|{}|{}\n", mv.getM2().maxCoeff(), kPrecision, min_partition_size, Unlocked);
//         assert(PQputCopyData(conn, data.c_str(), (int) data.length()) == 1);
//       }
//       #pragma omp for nowait
//       for (long long i = group_count; i < total_groups; i ++){
//         string data = fmt::format("{}|", i+1);
//         for (int j = 0; j < m; j ++){
//           data += fmt::format("{}|{}|", kNullLiteral, kNullLiteral);
//         }
//         data += fmt::format("{}|{}|{}\n", kNullLiteral, min_partition_size, Unitialized);
//         assert(PQputCopyData(conn, data.c_str(), (int) data.length()) == 1);
//       }
//       assert(PQputCopyEnd(conn, NULL) == 1);
//       res = PQgetResult(conn);
//       assert(PQresultStatus(res) == PGRES_COMMAND_OK);
//       PQclear(res);

//       // Populate P tabless
//       #pragma omp for
//       for (int i = 0; i < group_count; i ++){
//         const auto& g = *groups[i];
//         for (int j : g) group_indices[j] = i;
//       }
//       sql = fmt::format("COPY \"{}\" FROM STDIN with(delimiter ',');", p_name);
//       res = PQexec(conn, sql.c_str());
//       assert(PQresultStatus(res) == PGRES_COPY_IN);
//       PQclear(res);
//       int seg = omp_get_thread_num();
//       long long start_count = seg * overall_chunk + 1;
//       long long end_count = min((seg+1) * overall_chunk, n);
//       long long in_memory_count = ceilDiv(end_count-start_count+1, kInMemorySize);
//       for (long long i = 0; i < in_memory_count; i ++){
//         long long left = start_count + i*kInMemorySize;
//         long long right = min(left + kInMemorySize - 1, end_count);
//         string data = "";
//         for (long long tid = left; tid <= right; tid ++){
//           long long gid = 0;
//           if (tid <= kInitialSize) gid = group_indices[tid-1] + 1;
//           data += fmt::format("{},{}\n", tid, gid);
//         }
//         assert(PQputCopyData(conn, data.c_str(), (int) data.length()) == 1);
//       }
//       assert(PQputCopyEnd(conn, NULL) == 1);
//       res = PQgetResult(conn);
//       assert(PQresultStatus(res) == PGRES_COMMAND_OK);
//       PQclear(res);
//       PQfinish(conn);

//       #pragma omp barrier
//       #pragma omp for nowait
//       for (int i = 0; i < group_count; i ++) delete groups[i];
//     }
//     #pragma omp barrier

//     #pragma omp master
//     {
//       pro.stop(5);
//       pro.clock(6);
//       delete[] A;

//       string conninfo = fmt::format("postgresql://{}@{}?port={}&dbname={}&password={}", kPgUser, kPgHostaddr, kPgPort, dbname, kPgPassword);
//       PGconn *conn = PQconnectdb(conninfo.c_str());
//       assert(PQstatus(conn) == CONNECTION_OK);
//       PGresult *res = NULL;

//       pro.clock(8);
//       vector<string> intervals;
//       for (int i = 0; i < min(m, kMaxMultiColumnIndexes); i ++) intervals.push_back("interval_" + cols[i]);
//       string sql = fmt::format("CREATE INDEX \"{}_group_interval\" ON \"{}\" USING gist ({});", g_name, g_name, join(intervals, ","));
//       res = PQexec(conn, sql.c_str());
//       PQclear(res);
//       pro.stop(8);
//       pro.clock(9);
//       sql = fmt::format("ALTER TABLE \"{}\" ADD PRIMARY KEY ({});", g_name, kId);
//       res = PQexec(conn, sql.c_str());
//       PQclear(res);
//       pro.stop(9);
//       // pro.clock(10);
//       // sql = fmt::format("CREATE INDEX \"{}_tid_index\" ON \"{}\"(tid)", p_name, p_name);
//       // res = PQexec(conn, sql.c_str());
//       // PQclear(res);
//       // pro.stop(10);
//       pro.clock(11);
//       sql = fmt::format("CREATE INDEX \"{}_gid_index\" ON \"{}\"(gid) USING btree", p_name, p_name);
//       res = PQexec(conn, sql.c_str());
//       PQclear(res);
//       pro.stop(11);
//       // sql = fmt::format(""\
//       //   "CREATE OR REPLACE FUNCTION traverse( "\
//       //   "	start_id BIGINT, "\
//       //   "	end_id BIGINT "\
//       //   ") "\
//       //   "RETURNS TABLE( "\
//       //   "	v BIGINT "\
//       //   ") "\
//       //   "LANGUAGE plpgsql AS "\
//       //   "$$ "\
//       //   "DECLARE "\
//       //   "	query text := 'SELECT {},{} FROM \"{}\" WHERE {} BETWEEN %L AND %L'; "\
//       //   "	subquery text := 'SELECT {},max_size,status FROM \"{}\" AS t WHERE {} LIMIT 1'; "\
//       //   "	p_update_query text := 'UPDATE \"{}\" SET gid = %L WHERE tid = %L'; "\
//       //   "	p_size_query text := 'SELECT COUNT(*) FROM \"{}\" WHERE gid = %L'; "\
//       //   "	g_status_query text := 'SELECT status FROM \"{}\" WHERE {} = %L'; "\
//       //   "	g_lock_query text := 'UPDATE \"{}\" SET status = {}, max_size = max_size*2 WHERE {} = %L'; "\
//       //   "	subrec record; "\
//       //   "	rec record; "\
//       //   "	subrec_size BIGINT; "\
//       //   " subrec_status INTEGER; "\
//       //   "	is_locked bool; "\
//       //   "BEGIN "\
//       //   "	FOR rec IN EXECUTE format(query, start_id, end_id) LOOP "\
//       //   "		FOR subrec in EXECUTE format(subquery, {}) LOOP "\
//       //   "			IF	subrec.status = {} THEN "\
//       //   "				EXECUTE format(p_update_query, subrec.{}, rec.{}); "\
//       //   "				EXECUTE format(p_size_query, subrec.{}) INTO subrec_size; "\
//       //   "				IF subrec_size >= subrec.max_size THEN "\
//       //   "					SELECT pg_try_advisory_lock(subrec.{}) INTO is_locked; "\
//       //   "					EXECUTE format(g_status_query, subrec.{}) INTO subrec_status; "\
//       //   "					IF is_locked AND subrec_status = {} THEN "\
//       //   "						EXECUTE format(g_lock_query, subrec.{}); "\
//       //   "						v := rec.{}; "\
//       //   "						RETURN NEXT; "\
//       //   "						v := -subrec.{}; "\
//       //   "						RETURN NEXT; "\
//       //   "						PERFORM pg_advisory_unlock(subrec.{}); "\
//       //   "						RETURN; "\
//       //   "					ELSE "\
//       //   "						EXECUTE format(p_update_query, 0, rec.{}); "\
//       //   "						v := rec.{}; "\
//       //   "						RETURN NEXT; "\
//       //   "					END IF; "\
//       //   "				END IF; "\
//       //   "			ELSE "\
//       //   "				v := -rec.{}; "\
//       //   "				RETURN NEXT; "\
//       //   "			END IF; "\
//       //   "		END LOOP; "\
//       //   "	END LOOP; "\
//       //   "END; "\
//       //   "$$; ", kId, col_names, table_name, kId, kId, g_name, condition_names, p_name, p_name, g_name, kId, g_name, Locked, kId, rec_names, Unlocked, kId, kId, kId, kId, kId, Unlocked, kId, kId, kId, kId, kId, kId, kId);
//       // res = PQexec(conn, sql.c_str());
//       // PQclear(res);

//       sql = fmt::format("SELECT {}();", kClearLockFunction);
//       res = PQexec(conn, sql.c_str());
//       PQclear(res);
//       pro.stop(6);
//       PQfinish(conn);
//     }
//     #pragma omp barrier

//   //   // Phase 2
//   //   {
//   //     string conninfo = fmt::format("postgresql://{}@{}?port={}&dbname={}&password={}", kPgUser, kPgHostaddr, kPgPort, dbname, kPgPassword);
//   //     PGconn *conn = PQconnectdb(conninfo.c_str());
//   //     assert(PQstatus(conn) == CONNECTION_OK);
//   //     PGresult *res = NULL;
//   //     int seg = omp_get_thread_num();
//   //     long long start_id = kInitialSize + seg * phase_two_chunk + 1;
//   //     long long end_id = min(kInitialSize + (seg + 1) * phase_two_chunk, n);
//   //     vector<long long> remained_indices (kReserveSize);
//   //     int remained_size = 0;
//   //     string sql = fmt::format("SELECT * FROM {}({}, {});", kTraverseFunction, start_id, end_id);
//   //     #pragma omp critical
//   //     res = PQexec(conn, sql.c_str());
//   //     int ntuples = PQntuples(res);
//   //     long long group_id = 0;
//   //     long long last_id = 0;
//   //     if (ntuples >= 2){
//   //       group_id = atoll(PQgetvalue(res, ntuples-1, 0));
//   //       if (group_id < 0){
//   //         last_id = atoll(PQgetvalue(res, ntuples-2, 0));
//   //         ntuples -= 2;
//   //       }
//   //     }
//   //     for (int i = 0; i < ntuples; i ++){
//   //       remained_indices[remained_size] = atoll(PQgetvalue(res, i, 0));
//   //       remained_size ++;
//   //     }
//   //     PQclear(res);

//   //     if (remained_size >= kReserveSize){
//   //       vector<string> indices (remained_size);
//   //       for (int i = 0; i < remained_size; i ++) indices[i] = to_string(remained_indices[i]);
//   //       string index_names = join(indices, ",");
//   //       sql = fmt::format(""\
//   //         "CREATE OR REPLACE FUNCTION {}{}() "\
//   //         "RETURNS TABLE( "\
//   //         "	v BIGINT "\
//   //         ") "\
//   //         "LANGUAGE plpgsql AS "\
//   //         "$$ "\
//   //         "DECLARE "\
//   //         "	subquery text := 'SELECT {},status FROM \"{}\" AS t WHERE {} LIMIT 1'; "\
//   //         "	p_update_query text := 'UPDATE \"{}\" SET gid = %L WHERE tid = %L'; "\
//   //         "	subrec record; "\
//   //         "	rec record; "\
//   //         "BEGIN "\
//   //         "	FOR rec IN SELECT {},{} FROM \"{}\" WHERE {} IN ({}) LOOP "\
//   //         "		FOR subrec in EXECUTE format(subquery, {}) LOOP "\
//   //         "			IF subrec.status = {} THEN "\
//   //         "				EXECUTE format(p_update_query, subrec.{}, rec.{}); "\
//   //         "			ELSIF subrec.status = {} THEN "\
//   //         "				v := rec.{}; "\
//   //         "				RETURN NEXT; "\
//   //         "			END IF; "\
//   //         "		END LOOP; "\
//   //         "	END LOOP; "\
//   //         "END; "\
//   //         "$$; ", kProcessReserveFunction, seg, kId, g_name, condition_names, p_name, kId, col_names, table_name, kId, index_names, rec_names, Unlocked, kId, kId, Locked, kId);
//   //       res = PQexec(conn, sql.c_str());
//   //       PQclear(res);

//   //       sql = fmt::format("SELECT * FROM {}{}();", kProcessReserveFunction, seg);
//   //       res = PQexec(conn, sql.c_str());
//   //       remained_size = PQntuples(res);
//   //       for (int i = 0; i < remained_size; i ++) remained_indices[i] = atoll(PQgetvalue(res, i, 0));
//   //       PQclear(res);
//   //     }

//   //     if (group_id < 0){
//   //       sql = fmt::format("SELECT t.{},{} FROM \"{}\" AS t "\
//   //         "INNER JOIN \"{}\" AS p ON (t.{}=p.tid) WHERE p.gid={};", kId, t_names, table_name, p_name, kId, -group_id);
//   //       res = PQexec(conn, sql.c_str());
//   //       int size = PQntuples(res);
//   //       vector<long long> tids (size);
//   //       double *A = new double [size * m];
//   //       MeanVar mv = MeanVar(m);
//   //       for (int i = 0; i < size; i ++){
//   //         tids[i] = atoll(PQgetvalue(res, i, 0));
//   //         for (int j = 0; j < m; j ++) A[at(j, i)] = atof(PQgetvalue(res, i, j+1));
//   //         mv.add(A+at(0, i));
//   //       }
//   //       PQclear(res);
//   //       VectorXd total_vars = mv.getM2();
//   //       int mi = -1;
//   //       double max_total_var = -1;
//   //       for (int i = 0; i < m; i ++){
//   //         double total_var = total_vars(i);
//   //         if (max_total_var < total_var){
//   //           max_total_var = total_var;
//   //           mi = i;
//   //         }
//   //       }
//   //       double reduced_var = max_total_var / size * var_ratio;
//   //       #pragma omp critical
//   //       {
//   //         for (int i = 0; i < remained_size; i ++) cout << remained_indices[i] << " ";
//   //         cout << endl;
//   //         cout << group_count << " " << group_id << " " << size << " " << start_id << " " << last_id << endl;
//   //       }
//   //       // VectorXi indices (size); 
//   //       // iota(indices.begin(), indices.end(), 0);
//   //       // print(indices);
//   //       delete[] A;
//   //     }
//   //     PQfinish(conn);
//   //   }
//   }
// }