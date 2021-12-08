#define FMT_HEADER_ONLY

#include "fmt/core.h"
#include "dlv.h"
#include "libpq-fe.h"
#include "pb/util/uconfig.h"
#include "pb/util/umisc.h"

double DynamicLowVariance::kGroupRatio = 0.01;
int DynamicLowVariance::kMinPartitionSize = 1000000;
int DynamicLowVariance::kInitialSize = 100000;

void DynamicLowVariance::init(){
  pro.clock(0);
  string sql;
  string conninfo = fmt::format("postgresql://{}@{}?port={}&dbname={}&password={}", kPgUser, kPgHostaddr, kPgPort, dbname, kPgPassword);
  PGconn* conn = PQconnectdb(conninfo.c_str());
  assert(PQstatus(conn) == CONNECTION_OK);
  PGresult *res = NULL;
  PGresult *nres = NULL;

  sql = fmt::format("SET client_min_messages = warning;");
  res = PQexec(conn, sql.c_str());
  PQclear(res);

  sql = fmt::format(""
    "CREATE TABLE IF NOT EXISTS {}("
    "	table_name VARCHAR(63) NOT NULL,"
    "	count BIGINT DEFAULT 0,"
    " column_names text[],"
    "	mean DOUBLE PRECISION[],"
    "	M2 DOUBLE PRECISION[]"
    ")", kStatTable);
  res = PQexec(conn, sql.c_str());
  PQclear(res);

  sql = fmt::format(""
    "CREATE TABLE IF NOT EXISTS {}("
    "	table_name VARCHAR(63) NOT NULL,"
    "	partition_name VARCHAR(31) NOT NULL,"
    "	column_names text[],"
    "	group_ratio DOUBLE PRECISION,"
    "	min_partition_size INTEGER,"
    "	initial_size INTEGER,"
    " layer_count INTEGER"
    ");", kPartitionTable);
  res = PQexec(conn, sql.c_str());
  PQclear(res);

  string table_name = "test_data";
  sql = fmt::format("SELECT * FROM {} LIMIT 1;", table_name);
  res = PQexec(conn, sql.c_str());
  int col_count = PQnfields(res);
  vector<string> numeric_cols;
  for (int i = 0; i < col_count; i ++){
    string col_name (PQfname(res, i));
    if (!col_name.compare(kId)) continue;
    int oid = PQftype(res, i);
    sql = fmt::format("SELECT typcategory FROM pg_type WHERE oid={};", oid);
    nres = PQexec(conn, sql.c_str());
    char typ = PQgetvalue(nres, 0, 0)[0];
    PQclear(nres);
    if (typ == 'N') numeric_cols.push_back(col_name);
  }
  PQclear(res);

  sql = fmt::format("DROP FUNCTION IF EXISTS compute_stat;");
  res = PQexec(conn, sql.c_str());
  PQclear(res);

  string inject = ""
    "		delta := rec.{} - means[{}];"
    "		means[{}] := means[{}] + delta / count;"
    "		M2s[{}] := M2s[{}] + delta * (rec.{} - means[{}]);"
    "";
  string injects = "";
  for (int i = 1; i <= numeric_cols.size(); i ++){
    string col = numeric_cols[i-1];
    injects += fmt::format(inject, col, i, i, i, i, i, col, i);
  }

  sql = fmt::format(""
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
    "$$;", join(numeric_cols, ","), table_name, kId, numeric_cols.size(), numeric_cols.size(), injects, numeric_cols.size());
  res = PQexec(conn, sql.c_str());
  PQclear(res);

  sql = fmt::format("SELECT {} FROM {} ORDER BY {} DESC LIMIT 1;", kId, table_name, kId);
  res = PQexec(conn, sql.c_str());
  int count = atoi(PQgetvalue(res, 0, 0)) + 1;
  cout << count << '\n';
  PQclear(res);
  // #pragma omp parallel num_threads(kPCore)
  // {

  // }
  int start_id = 0;
  int end_id = 124999;
  sql = fmt::format("SELECT * FROM compute_stat({}, {});", start_id, end_id);
  res = PQexec(conn, sql.c_str());
  VectorXd mean (numeric_cols.size());
  VectorXd M2 (numeric_cols.size());
  for (int i = 0; i < numeric_cols.size(); i ++){
    mean(i) = atof(PQgetvalue(res, i, 0));
    M2(i) = atof(PQgetvalue(res, i, 1));
  }
  print(mean);
  print(M2);
  PQclear(res);

  pro.stop(0);
}

DynamicLowVariance::DynamicLowVariance(string dbname): dbname(dbname){
  vector<string> names = {"Init", "1", "2", "3", "4", "5", "6"};
  pro = Profiler(names);
  init();
}