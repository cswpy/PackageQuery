#include "det_sql.h"

#include "pb/util/uconfig.h"
#include "pb/util/unumeric.h"
#include "pb/util/umisc.h"
#include "pb/lib/random_quantile.h"

#include "libpq-fe.h"

const double kMaxLog2OfInt = 31.0;

DetSql::~DetSql(){
  delete pg;
}

DetSql::DetSql(string table_name, string obj_col, bool is_maximize, vector<string> att_cols, vector<int> att_senses, bool has_count_constraint, long long u)
  : table_name(table_name), obj_col(obj_col), is_maximize(is_maximize), att_cols(att_cols), att_senses(att_senses), u(u), has_count_constraint(has_count_constraint){
  pg = new PgManager();
  table_size = pg->getSize(table_name);
  table_cols = pg->listColumns(table_name);
  assert(isIn(table_cols, obj_col));
  for (string col : att_cols) assert(isIn(table_cols, col));
}

void DetSql::addFilter(string col, double l, double u){
  if (!isIn(table_cols, col)) return;
  if (l == -DBL_MAX && u == DBL_MAX) return;
  filter_cols.push_back(col);
  filter_intervals.emplace_back(l, u);
}

void DetSql::addFilterWithRatio(string col, double F, int sense){
  long long size = pg->getSize(table_name);
  int core = kPCore;
  double main_memory = kMainMemorySize;
  long long chunk = ceilDiv(size, (long long) core);
  long long slide_size = getTupleCount(1, core+1, main_memory);
  double eps = kMaxLog2OfInt/slide_size;
  RandomQuantile rq (eps);
  #pragma omp parallel num_threads(core)
  {
    string sql;
    PGconn *conn = PQconnectdb(pg->conninfo.c_str());
    assert(PQstatus(conn) == CONNECTION_OK);
    PGresult *res = NULL;

    int seg = omp_get_thread_num();
    long long start_id = seg * chunk + 1;
    long long end_id = min((seg + 1) * chunk, size);
    long long slide_count = ceilDiv(end_id - start_id + 1, slide_size);
    for (long long sl = 0; sl < slide_count; sl ++){
      long long slide_start_id = start_id + sl * slide_size;
      long long slide_end_id = min(start_id + (sl + 1) * slide_size - 1, end_id);
      sql = fmt::format("SELECT {} FROM \"{}\" WHERE {} BETWEEN {} AND {};", col, table_name, kId, slide_start_id, slide_end_id);
      res = PQexec(conn, sql.c_str());
      vector<double> vals (PQntuples(res));
      for (long long i = 0; i < PQntuples(res); i++){
        vals[i] = atof(PQgetvalue(res, i, 0));
      }
      PQclear(res);
      #pragma omp critical (c0)
      {
        for (auto &v : vals) rq.feed(v);
      }
    }
    PQfinish(conn);
  }
  rq.finalize();
  switch(sense) {
    case LowerBounded:
      addFilter(col, rq.query_for_value(1-F), DBL_MAX);
      break;
    case UpperBounded:
      addFilter(col, -DBL_MAX, rq.query_for_value(F));
      break;
    case Bounded:
    ///// (1-F)/2 //// l //////F////// u ///// (1-F)/2 //////
      addFilter(col, rq.query_for_value((1-F)/2), rq.query_for_value(1-(1-F)/2));
      break;
    default:
      throw invalid_argument("Invalid DetSense type: " + to_string(sense));
  }
}
  
bool DetSql::isFiltering(){
  return filter_cols.size() > 0;
}