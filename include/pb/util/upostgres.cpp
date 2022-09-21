#define FMT_HEADER_ONLY

#include <iostream>
#include <omp.h>
#include <cfloat>
#include "fmt/core.h"
#include "upostgres.h"
#include "pb/util/uconfig.h"
#include "pb/util/umisc.h"
#include "pb/util/unumeric.h"

using std::min;
using std::max;
using std::cout;
using std::endl;

// const string kComputeStatFunction = "compute_stat";

void showError(PGconn *conn){
  string error (PQerrorMessage(conn));
  if (error.length()) cout << error << endl;
  else cout << "No errors" << endl;
}

PgManager::~PgManager(){
  PQfinish(_conn);
}

PgManager::PgManager(){
  conninfo = fmt::format("postgresql://{}@{}?port={}&dbname={}&password={}", kPgUser, kPgHostaddr, kPgPort, kPgDatabase, kPgPassword);
  _conn = PQconnectdb(conninfo.c_str());
  assert(PQstatus(_conn) == CONNECTION_OK);
  _res = NULL;

  _sql = fmt::format("SET client_min_messages = warning;");
  _res = PQexec(_conn, _sql.c_str());
  PQclear(_res);
}

long long PgManager::getSize(string table_name){
  _sql = fmt::format("SELECT {} FROM \"{}\" ORDER BY {} DESC LIMIT 1;", kId, table_name, kId);
  _res = PQexec(_conn, _sql.c_str());
  long long size = atoll(PQgetvalue(_res, 0, 0));
  PQclear(_res);
  return size;
}

vector<string> PgManager::getNumericCols(string table_name){
  vector<string> cols;
  _sql = fmt::format("SELECT * FROM \"{}\" LIMIT 1;", table_name);
  _res = PQexec(_conn, _sql.c_str());
  int col_count = PQnfields(_res);
  PGresult *res = NULL;
  for (int i = 0; i < col_count; i ++){
    string col_name (PQfname(_res, i));
    if (!col_name.compare(kId)) continue;
    int oid = PQftype(_res, i);
    _sql = fmt::format("SELECT typcategory FROM pg_type WHERE oid={};", oid);
    res = PQexec(_conn, _sql.c_str());
    char typ = PQgetvalue(res, 0, 0)[0];
    PQclear(res);
    if (typ == 'N') cols.push_back(col_name);
  }
  PQclear(_res);
  return cols;
}

Stat* PgManager::computeStats(string table_name, const vector<string> &cols){
  long long size = getSize(table_name);
  long long chunk = ceilDiv(size, kPCore);
  Stat* stat = new Stat(cols);
  #pragma omp parallel num_threads(kPCore)
  {
    string sql;
    PGconn* conn = PQconnectdb(conninfo.c_str());
    assert(PQstatus(conn) == CONNECTION_OK);
    PGresult *res = NULL;

    int seg = omp_get_thread_num();
    long long start_id = seg * chunk + 1;
    long long end_id = min((seg + 1) * chunk, size);
    MeanVar mv = MeanVar(cols.size());

    sql = fmt::format("SELECT {} FROM \"{}\" WHERE {} BETWEEN {} AND {};", join(cols, ","), table_name, kId, start_id, end_id);
    res = PQexec(conn, sql.c_str());
    VectorXd x (cols.size());
    VectorXd amin (cols.size()); amin.fill(DBL_MAX);
    VectorXd amax (cols.size()); amax.fill(-DBL_MAX);
    for (int i = 0; i < PQntuples(res); i++){
      for (int j = 0; j < (int) cols.size(); j ++){
        x(j) = atof(PQgetvalue(res, i, j));
        amin(j) = min(amin(j), x(j));
        amax(j) = max(amax(j), x(j));
      }
      mv.add(x);
    }

    #pragma omp critical
    {
      stat->add(mv.sample_count, mv.getMean(), mv.getM2(), amin, amax);
    }
  }
  return stat;
}

bool PgManager::existTable(string table_name){
  _sql = fmt::format("SELECT * FROM pg_tables WHERE tablename='{}' AND schemaname='public'", table_name);
  _res = PQexec(_conn, _sql.c_str());
  bool exists = PQntuples(_res) > 0;
  PQclear(_res);
  return exists;
}

void PgManager::dropTable(string table_name){
  _sql = fmt::format("DROP TABLE IF EXISTS \"{}\";", table_name);
  _res = PQexec(_conn, _sql.c_str());
  PQclear(_res);
}

// Deprecated since Postgres does not allow inter-multithreaded

// Stat* PgManager::computeStats(string table_name, const vector<string> &cols){
//   string inject = ""
//     "		delta := rec.{} - means[{}];"
//     "		means[{}] := means[{}] + delta / count;"
//     "		M2s[{}] := M2s[{}] + delta * (rec.{} - means[{}]);"
//     "   amins[{}] := LEAST(amins[{}], rec.{});"
//     "   amaxs[{}] := GREATEST(amaxs[{}], rec.{});"
//     "";
//   string injects = "";
//   for (int i = 1; i <= (int) cols.size(); i ++){
//     string col = cols[i-1];
//     injects += fmt::format(inject, col, i, i, i, i, i, col, i, i, i, col, i, i, col);
//   }

//   _sql = fmt::format("DROP FUNCTION {};", kComputeStatFunction);
//   _res = PQexec(_conn, _sql.c_str());
//   PQclear(_res);

//   _sql = fmt::format(""
//     "CREATE FUNCTION {}("
//     "	start_id BIGINT,"
//     "	end_id BIGINT"
//     ")"
//     "RETURNS TABLE("
//     "	mean DOUBLE PRECISION,"
//     "	M2 DOUBLE PRECISION,"
//     "	amin DOUBLE PRECISION,"
//     "	amax DOUBLE PRECISION"
//     ")"
//     "LANGUAGE plpgsql AS "
//     "$$"
//     "DECLARE"
//     "	rec record;"
//     "	query text := 'SELECT {} FROM \"{}\" WHERE {} BETWEEN %L AND %L;';"
//     "	count BIGINT := 0;"
//     "	means double precision[] := array_fill(0, ARRAY[{}]);"
//     "	M2s double precision[] := array_fill(0, ARRAY[{}]);"
//     "	amins double precision[] := array_fill(1e300, ARRAY[{}]);"
//     "	amaxs double precision[] := array_fill(-1e300, ARRAY[{}]);"
//     "	delta double precision;"
//     "BEGIN"
//     "	FOR rec IN EXECUTE format(query, start_id, end_id) LOOP"
//     "		count := count + 1;"
//     "{}"
//     "	END LOOP;"
//     "	FOR i in 1..{} LOOP"
//     "		mean := means[i];"
//     "		M2 := M2s[i];"
//     "		amin := amins[i];"
//     "		amax := amaxs[i];"
//     "		RETURN NEXT;"
//     "	END LOOP;"
//     "END;"
//     "$$;", kComputeStatFunction, join(cols, ","), table_name, kId, cols.size(), cols.size(), cols.size(), cols.size(), injects, cols.size());
//   _res = PQexec(_conn, _sql.c_str());
//   PQclear(_res);
//   long long size = getSize(table_name);
//   long long chunk = ceilDiv(size, kPCore);
//   Stat* stat = new Stat(cols);
//   #pragma omp parallel num_threads(kPCore)
//   {
//     string sql;
//     string conninfo = fmt::format("postgresql://{}@{}?port={}&dbname={}&password={}", kPgUser, kPgHostaddr, kPgPort, dbname, kPgPassword);
//     PGconn* conn = PQconnectdb(conninfo.c_str());
//     assert(PQstatus(conn) == CONNECTION_OK);
//     PGresult *res = NULL;

//     int seg = omp_get_thread_num();
//     long long start_id = seg * chunk + 1;
//     long long end_id = min((seg + 1) * chunk, size);
//     sql = fmt::format("SELECT * FROM {}({}, {});", kComputeStatFunction, start_id, end_id);
//     res = PQexec(conn, sql.c_str());
//     #pragma omp critical
//     {
//       cout << start_id << " " << end_id << endl;
//     }
//     VectorXd local_mean (cols.size());
//     VectorXd local_M2 (cols.size());
//     VectorXd local_amin (cols.size());
//     VectorXd local_amax (cols.size());
//     for (int i = 0; i < (int) cols.size(); i ++){
//       local_mean(i) = atof(PQgetvalue(res, i, 0));
//       local_M2(i) = atof(PQgetvalue(res, i, 1));
//       local_amin(i) = atof(PQgetvalue(res, i, 2));
//       local_amax(i) = atof(PQgetvalue(res, i, 3));
//     }
//     PQclear(res);
//     PQfinish(conn);
//     #pragma omp critical
//     {
//       stat->add(end_id-start_id+1, local_mean, local_M2, local_amin, local_amax);
//     }
//   }
//   return stat;
// }

Stat::Stat(vector<string> cols): cols(cols){
  size = 0LL;
  mean.resize(cols.size()); mean.fill(0);
  M2.resize(cols.size()); M2.fill(0);
  amin.resize(cols.size()); amin.fill(DBL_MAX);
  amax.resize(cols.size()); amax.fill(-DBL_MAX);
}

void Stat::add(long long size, VectorXd mean, VectorXd M2){
  if (!this->size){
    this->size = size;
    this->mean = mean;
    this->M2 = M2;
  } else{
    long long tmp_size = this->size + size; 
    VectorXd delta = mean - this->mean;
    this->mean = (this->size * this->mean + size * mean) / tmp_size;
    this->M2 += M2 + delta.cwiseProduct(delta) * (this->size * size / (long double) (tmp_size));
    this->size = tmp_size;
  }
}

void Stat::add(long long size, VectorXd mean, VectorXd M2, VectorXd amin, VectorXd amax){
  for (int i = 0; i < (int) cols.size(); i ++){
    this->amin(i) = min(this->amin(i), amin(i));
    this->amax(i) = max(this->amax(i), amax(i));
  }
  add(size, mean, M2);
}

int Stat::getIndex(string col){
  for (int i = 0; i < (int) cols.size(); i ++){
    if (!cols[i].compare(col)) return i;
  }
  return -1;
}

double Stat::getVar(int i){
  assert(0 <= i && i < (int) cols.size());
  if (size <= 1) return M2(i);
  return M2(i) / (size - 1);
}