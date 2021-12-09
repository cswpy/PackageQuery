#define FMT_HEADER_ONLY

#include "fmt/core.h"
#include "upostgres.h"
#include "pb/util/uconfig.h"
#include "pb/util/umisc.h"
#include "pb/util/unumeric.h"

using std::min;

PgManager::~PgManager(){
  PQfinish(_conn);
}

PgManager::PgManager(string dbname): dbname(dbname){
  string conninfo = fmt::format("postgresql://{}@{}?port={}&dbname={}&password={}", kPgUser, kPgHostaddr, kPgPort, dbname, kPgPassword);
  _conn = PQconnectdb(conninfo.c_str());
  assert(PQstatus(_conn) == CONNECTION_OK);
  _res = NULL;

  _sql = fmt::format("SET client_min_messages = warning;");
  _res = PQexec(_conn, _sql.c_str());
  PQclear(_res);
}

void PgManager::createExtention(string extension){
  _sql = fmt::format("CREATE EXTENSION IF NOT EXISTS {};", extension);
  _res = PQexec(_conn, _sql.c_str());
  PQclear(_res);
}

long long PgManager::getSize(string table_name){
  _sql = fmt::format("SELECT {} FROM {} ORDER BY {} DESC LIMIT 1;", kId, table_name, kId);
  _res = PQexec(_conn, _sql.c_str());
  long long size = atoll(PQgetvalue(_res, 0, 0));
  PQclear(_res);
  return size;
}

vector<string> PgManager::getNumericCols(string table_name){
  vector<string> cols;
  _sql = fmt::format("SELECT * FROM {} LIMIT 1;", table_name);
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

Stat* PgManager::computeStats(string table_name){
  return computeStats(table_name, getNumericCols(table_name));
}

Stat* PgManager::computeStats(string table_name, const vector<string> &cols){
  _sql = fmt::format("DROP FUNCTION IF EXISTS compute_stat;");
  _res = PQexec(_conn, _sql.c_str());
  PQclear(_res);

  string inject = ""
    "		delta := rec.{} - means[{}];"
    "		means[{}] := means[{}] + delta / count;"
    "		M2s[{}] := M2s[{}] + delta * (rec.{} - means[{}]);"
    "";
  string injects = "";
  for (int i = 1; i <= cols.size(); i ++){
    string col = cols[i-1];
    injects += fmt::format(inject, col, i, i, i, i, i, col, i);
  }

  _sql = fmt::format(""
    "CREATE OR REPLACE FUNCTION compute_stat("
    "	start_id BIGINT,"
    "	end_id BIGINT"
    ")"
    "RETURNS TABLE("
    "	mean DOUBLE PRECISION,"
    "	M2 DOUBLE PRECISION"
    ")"
    "LANGUAGE plpgsql AS "
    "$$"
    "DECLARE"
    "	rec record;"
    "	query text := 'SELECT {} FROM {} WHERE {} BETWEEN %L AND %L;';"
    "	count BIGINT := 0;"
    "	means double precision[] := array_fill(0, ARRAY[{}]);"
    "	M2s double precision[] := array_fill(0, ARRAY[{}]);"
    "	delta double precision;"
    "BEGIN"
    "	FOR rec IN EXECUTE format(query, start_id, end_id) LOOP"
    "		count := count + 1;"
    "{}"
    "	END LOOP;"
    "	FOR i in 1..{} LOOP"
    "		mean := means[i];"
    "		M2 := M2s[i];"
    "		RETURN NEXT;"
    "	END LOOP;"
    "END;"
    "$$;", join(cols, ","), table_name, kId, cols.size(), cols.size(), injects, cols.size());
  _res = PQexec(_conn, _sql.c_str());
  PQclear(_res);

  long long size = getSize(table_name);
  long long chunk = ceilDiv(size, kPCore);
  Stat* stat = new Stat(cols);
  #pragma omp parallel num_threads(kPCore)
  {
    string sql;
    string conninfo = fmt::format("postgresql://{}@{}?port={}&dbname={}&password={}", kPgUser, kPgHostaddr, kPgPort, dbname, kPgPassword);
    PGconn* conn = PQconnectdb(conninfo.c_str());
    assert(PQstatus(conn) == CONNECTION_OK);
    PGresult *res = NULL;

    int seg = omp_get_thread_num();
    long long start_id = seg * chunk + 1;
    long long end_id = min((seg + 1) * chunk, size);
    sql = fmt::format("SELECT * FROM compute_stat({}, {});", start_id, end_id);
    res = PQexec(conn, sql.c_str());
    VectorXd local_mean (cols.size());
    VectorXd local_M2 (cols.size());
    for (int i = 0; i < cols.size(); i ++){
      local_mean(i) = atof(PQgetvalue(res, i, 0));
      local_M2(i) = atof(PQgetvalue(res, i, 1));
    }
    PQclear(res);
    PQfinish(conn);
    long long local_count = end_id - start_id + 1;
    #pragma omp critical
    {
      stat->add(end_id-start_id+1, local_mean, local_M2);
    }
  }
  return stat;
}

Stat::Stat(vector<string> cols): cols(cols){
  size = 0LL;
  mean.resize(cols.size()); mean.fill(0);
  M2.resize(cols.size()); M2.fill(0);
}

void Stat::add(long long size, VectorXd &mean, VectorXd &M2){
  long long tmp_size = this->size + size; 
  VectorXd delta = mean - this->mean;
  this->mean = (this->size * this->mean + size * mean) / tmp_size;
  this->M2 += M2 + delta.cwiseProduct(delta) * (this->size * size / (long double) (tmp_size));
  this->size = tmp_size;
}