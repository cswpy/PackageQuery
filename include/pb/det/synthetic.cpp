#define FMT_HEADER_ONLY

#include "fmt/core.h"
#include "synthetic.h"
#include "libpq-fe.h"
#include "pb/util/uconfig.h"
#include "pb/util/unumeric.h"

using std::min;

string Synthetic::table_name = "test_data";
int Synthetic::in_memory_chunk = 1000000;

Synthetic::Synthetic(string dbname): dbname(dbname){
  vector<string> names = {"Init", "Create table", "Generate data", "Create index"};
  pro = Profiler(names);
  pro.clock(0);
  init();
  pro.stop(0);
}

void Synthetic::createUniform(long long N, double mean, double var){
  create("uniform_4col", N, {mean-sqrt(3*var), mean+sqrt(3*var)});
}

void Synthetic::createNormal(long long N, double mean, double var){
  create("normal_4col", N, {mean, sqrt(var)});
}

void Synthetic::create(string method, long long N, vector<double> args){
  long long chunk = ceilDiv(N, kPCore);
  pro.clock(1);
  {
    string sql;
    string conninfo = fmt::format("postgresql://{}@{}?port={}&dbname={}&password={}", kPgUser, kPgHostaddr, kPgPort, dbname, kPgPassword);
    PGconn* conn = PQconnectdb(conninfo.c_str());
    assert(PQstatus(conn) == CONNECTION_OK);
    PGresult *res = NULL;

    sql = fmt::format("DROP TABLE IF EXISTS {} CASCADE;", Synthetic::table_name);
    res = PQexec(conn, sql.c_str());
    PQclear(res);

    sql = fmt::format(""
      "CREATE TABLE IF NOT EXISTS {} ("
      "	id BIGINT,"
      "	a1 DOUBLE PRECISION,"
      "	a2 DOUBLE PRECISION,"
      "	a3 DOUBLE PRECISION,"
      "	a4 DOUBLE PRECISION"
      ");", Synthetic::table_name);
    res = PQexec(conn, sql.c_str());
    PQclear(res);

    PQfinish(conn);
  }
  pro.stop(1);
  pro.clock(2);
  #pragma omp parallel num_threads(kPCore)
  {
    string sql;
    string conninfo = fmt::format("postgresql://{}@{}?port={}&dbname={}&password={}", kPgUser, kPgHostaddr, kPgPort, dbname, kPgPassword);
    PGconn* conn = PQconnectdb(conninfo.c_str());
    assert(PQstatus(conn) == CONNECTION_OK);
    PGresult *res = NULL;
    int seg = omp_get_thread_num();
    long long start_count = seg * chunk + 1;
    long long end_count = min((seg + 1) * chunk, N);
    int in_memory_count = (int) ceilDiv(end_count-start_count+1, Synthetic::in_memory_chunk);
    for (int i = 0; i < in_memory_count; i ++){
      long long left = start_count + i * Synthetic::in_memory_chunk;
      long long right = min(left + Synthetic::in_memory_chunk - 1, end_count);
      sql = fmt::format("INSERT INTO {} SELECT * FROM {}({}, {}, {}, {});", Synthetic::table_name, method, left, right, args[0], args[1]);
      res = PQexec(conn, sql.c_str());
      PQclear(res);
    }
    PQfinish(conn);
  }
  pro.stop(2);
  pro.clock(3);
  {
    string sql;
    string conninfo = fmt::format("postgresql://{}@{}?port={}&dbname={}&password={}", kPgUser, kPgHostaddr, kPgPort, dbname, kPgPassword);
    PGconn* conn = PQconnectdb(conninfo.c_str());
    assert(PQstatus(conn) == CONNECTION_OK);
    PGresult *res = NULL;

    sql = fmt::format("ALTER TABLE {} ADD PRIMARY KEY (id);", Synthetic::table_name);
    res = PQexec(conn, sql.c_str());
    PQclear(res);

    PQfinish(conn);
  }
  pro.stop(3);
}

void Synthetic::init(){
  string sql;
  string conninfo = fmt::format("postgresql://{}@{}?port={}&dbname={}&password={}", kPgUser, kPgHostaddr, kPgPort, dbname, kPgPassword);
  PGconn* conn = PQconnectdb(conninfo.c_str());
  assert(PQstatus(conn) == CONNECTION_OK);
  PGresult *res = NULL;

  sql = ""
    "CREATE OR REPLACE FUNCTION uniform_4col("
    "	start_count BIGINT DEFAULT 1,"
    "    end_count BIGINT DEFAULT 1,"
    "    min DOUBLE PRECISION DEFAULT 0.0,"
    "    max DOUBLE PRECISION DEFAULT 1.0"
    ") "
    "RETURNS TABLE ("
    "	id BIGINT,"
    "	a1 DOUBLE PRECISION,"
    "	a2 DOUBLE PRECISION,"
    "	a3 DOUBLE PRECISION,"
    "	a4 DOUBLE PRECISION"
    ")"
    "RETURNS NULL ON NULL INPUT AS "
    "$$"
    "BEGIN"
    "	RETURN QUERY SELECT "
    "		UNNEST(ARRAY(SELECT GENERATE_SERIES FROM GENERATE_SERIES(start_count, end_count))) AS id,"
    "		UNNEST(ARRAY(SELECT min + (max - min) * RANDOM() FROM GENERATE_SERIES(0, end_count-start_count))) AS a1,"
    "		UNNEST(ARRAY(SELECT min + (max - min) * RANDOM() FROM GENERATE_SERIES(0, end_count-start_count))) AS a2,"
    "		UNNEST(ARRAY(SELECT min + (max - min) * RANDOM() FROM GENERATE_SERIES(0, end_count-start_count))) AS a3,"
    "		UNNEST(ARRAY(SELECT min + (max - min) * RANDOM() FROM GENERATE_SERIES(0, end_count-start_count))) AS a4;"
    "END;"
    "$$ LANGUAGE plpgsql;";
  res = PQexec(conn, sql.c_str());
  PQclear(res);

  sql = ""
    "CREATE OR REPLACE FUNCTION normal_4col("
    "	start_count BIGINT DEFAULT 1,"
    "    end_count BIGINT DEFAULT 1,"
    "    mean DOUBLE PRECISION DEFAULT 0.0,"
    "    stddev DOUBLE PRECISION DEFAULT 1.0"
    ") "
    "RETURNS TABLE ("
    "	id BIGINT,"
    "	a1 DOUBLE PRECISION,"
    "	a2 DOUBLE PRECISION,"
    "	a3 DOUBLE PRECISION,"
    "	a4 DOUBLE PRECISION"
    ")"
    "RETURNS NULL ON NULL INPUT AS "
    "$$"
    "BEGIN"
    "	RETURN QUERY SELECT "
    "		UNNEST(ARRAY(SELECT GENERATE_SERIES FROM GENERATE_SERIES(start_count, end_count))) AS id,"
    "		UNNEST(ARRAY(SELECT * FROM NORMAL_RAND((end_count-start_count+1)::Integer, mean, stddev))) AS a1,"
    "		UNNEST(ARRAY(SELECT * FROM NORMAL_RAND((end_count-start_count+1)::Integer, mean, stddev))) AS a2,"
    "		UNNEST(ARRAY(SELECT * FROM NORMAL_RAND((end_count-start_count+1)::Integer, mean, stddev))) AS a3,"
    "		UNNEST(ARRAY(SELECT * FROM NORMAL_RAND((end_count-start_count+1)::Integer, mean, stddev))) AS a4;"
    "END;"
    "$$ LANGUAGE plpgsql;";
  res = PQexec(conn, sql.c_str());
  PQclear(res);

  sql = "SET random_page_cost=1.1;";
  res = PQexec(conn, sql.c_str());
  PQclear(res);

  PQfinish(conn);
}