#define FMT_HEADER_ONLY

#include "fmt/core.h"
#include "synthetic.h"
#include "libpq-fe.h"
#include "pb/util/uconfig.h"
#include "pb/util/unumeric.h"
#include "pb/util/udeclare.h"
#include "pb/util/umisc.h"

using namespace pb;

string Synthetic::table_name = "test_data";

// const string kGenerateFunc = "generate_table";

Synthetic::Synthetic(string dbname): dbname(dbname){
  vector<string> names = {"Init", "Create table", "Generate data", "Create index"};
  pro = Profiler(names);
  pro.clock(0);
  init();
  pro.stop(0);
}

void Synthetic::createUniform(long long N, int count, double mean, double var){
  vector<double> means (count, mean);  
  vector<double> vars (count, var);
  create(N, count, 0, means, vars);
}

void Synthetic::createNormal(long long N, int count, double mean, double var){
  vector<double> means (count, mean);  
  vector<double> vars (count, var);
  create(N, 0, count, means, vars);
}

void Synthetic::createMixed(long long N, int ucount, int ncount, double mean, double var){
  int count = ucount + ncount;
  vector<double> means (count, mean);  
  vector<double> vars (count, var);
  create(N, ucount, ncount, means, vars);
}

void Synthetic::create(long long N, int ucount, int ncount, vector<double> means, vector<double> vars){
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

    string inject1 = ", u{} DOUBLE PRECISION";
    string inject2 = ", n{} DOUBLE PRECISION";
    string declare_inject = "";
    for (int i = 1; i <= ucount; i ++) declare_inject += fmt::format(inject1, i);
    for (int i = 1; i <= ncount; i ++) declare_inject += fmt::format(inject2, i);

    sql = fmt::format(""
      "CREATE TABLE IF NOT EXISTS {} ("
      "	{} BIGINT"
      " {}"
      ");", Synthetic::table_name, kId, declare_inject);
    res = PQexec(conn, sql.c_str());
    PQclear(res);
    PQfinish(conn);
  }
  pro.stop(1);

  pro.clock(2);
  long long chunk = ceilDiv(N, kPCore);
  #pragma omp parallel num_threads(kPCore)
  {
    int seg = omp_get_thread_num();
    std::default_random_engine gen;
    gen.seed(seg);
    std::uniform_real_distribution<double> u_dist(-sqrt(3), sqrt(3));
    std::normal_distribution<double> n_dist(0.0, 1.0);
    long long start_id = seg * chunk + 1;
    long long end_id = min((seg + 1) * chunk, N);

    string conninfo = fmt::format("postgresql://{}@{}?port={}&dbname={}&password={}", kPgUser, kPgHostaddr, kPgPort, dbname, kPgPassword);
    PGconn *conn = PQconnectdb(conninfo.c_str());
    assert(PQstatus(conn) == CONNECTION_OK);
    PGresult *res = NULL;

    string sql = fmt::format("COPY \"{}\" FROM STDIN with(delimiter ',');", Synthetic::table_name);
    res = PQexec(conn, sql.c_str());
    assert(PQresultStatus(res) == PGRES_COPY_IN);
    PQclear(res);
    VectorXd x (ucount + ncount);
    for (long long id=start_id; id <= end_id; id ++){
      for (int i = 0; i < ucount; i ++) x(i) = u_dist(gen) * sqrt(vars[i]) + means[i];
      for (int i = ucount; i < ucount + ncount; i ++) x(i) = n_dist(gen) * sqrt(vars[i]) + means[i];
      string data = fmt::format("{},{}\n", id, join(x, kPrecision));
      assert(PQputCopyData(conn, data.c_str(), (int) data.length()) == 1);
    }
    assert(PQputCopyEnd(conn, NULL) == 1);
    res = PQgetResult(conn);
    assert(PQresultStatus(res) == PGRES_COMMAND_OK);
    PQclear(res);
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

    sql = fmt::format("ALTER TABLE {} ADD PRIMARY KEY ({});", Synthetic::table_name, kId);
    res = PQexec(conn, sql.c_str());
    PQclear(res);

    PQfinish(conn);
  }
  pro.stop(3);
}

// Deprecated 
// void Synthetic::create(long long N, int ucount, int ncount, vector<double> means, vector<double> vars){
//   long long chunk = ceilDiv(N, kPCore);
//   pro.clock(1);
//   {
//     string sql;
//     string conninfo = fmt::format("postgresql://{}@{}?port={}&dbname={}&password={}", kPgUser, kPgHostaddr, kPgPort, dbname, kPgPassword);
//     PGconn* conn = PQconnectdb(conninfo.c_str());
//     assert(PQstatus(conn) == CONNECTION_OK);
//     PGresult *res = NULL;

//     sql = fmt::format("SET client_min_messages = warning;");
//     res = PQexec(conn, sql.c_str());
//     PQclear(res);

//     sql = fmt::format("DROP FUNCTION IF EXISTS {};", kGenerateFunc);
//     res = PQexec(conn, sql.c_str());
//     PQclear(res);

//     string inject1 = ", u{} DOUBLE PRECISION";
//     string inject2 = ", n{} DOUBLE PRECISION";
//     string declare_inject = "";
//     for (int i = 1; i <= ucount; i ++) declare_inject += fmt::format(inject1, i);
//     for (int i = 1; i <= ncount; i ++) declare_inject += fmt::format(inject2, i);

//     inject1 = ", UNNEST(ARRAY(SELECT {:.{}Lf} + {:.{}Lf} * RANDOM() FROM GENERATE_SERIES(0, end_id-start_id))) AS u{}";
//     inject2 = ", UNNEST(ARRAY(SELECT * FROM NORMAL_RAND((end_id-start_id+1)::Integer, {:.{}Lf}, {:.{}Lf}))) AS n{}";
//     string query_inject = "";
//     for (int i = 0; i < ucount; i ++){
//       double tmp = sqrt(3*vars[i]);
//       query_inject += fmt::format(inject1, means[i]-tmp, kPrecision, 2*tmp, kPrecision, i+1);
//     }
//     for (int i = 0; i < ncount; i ++) query_inject += fmt::format(inject2, means[i+ucount], kPrecision, sqrt(vars[i+ucount]), kPrecision, i+1);
//     sql = fmt::format(""
//       "CREATE OR REPLACE FUNCTION {}( "
//       "	start_id BIGINT, "
//       "	end_id BIGINT "
//       ") "
//       "RETURNS TABLE( "
//       "{} BIGINT"
//       "{}"
//       ") "
//       "LANGUAGE plpgsql AS "
//       "$$ "
//       "BEGIN "
//       "	RETURN QUERY SELECT "
//       "	UNNEST(ARRAY(SELECT GENERATE_SERIES FROM GENERATE_SERIES(start_id, end_id))) AS {}"
//       " {};"
//       "END; "
//       "$$; ", kGenerateFunc, kId, declare_inject, kId, query_inject);
//     res = PQexec(conn, sql.c_str());
//     PQclear(res);

//     sql = fmt::format("DROP TABLE IF EXISTS {} CASCADE;", Synthetic::table_name);
//     res = PQexec(conn, sql.c_str());
//     PQclear(res);

//     sql = fmt::format(""
//       "CREATE TABLE IF NOT EXISTS {} ("
//       "	{} BIGINT"
//       " {}"
//       ");", Synthetic::table_name, kId, declare_inject);
//     res = PQexec(conn, sql.c_str());
//     PQclear(res);

//     PQfinish(conn);
//   }
//   pro.stop(1);
//   pro.clock(2);
//   #pragma omp parallel num_threads(kPCore)
//   {
//     string sql;
//     string conninfo = fmt::format("postgresql://{}@{}?port={}&dbname={}&password={}", kPgUser, kPgHostaddr, kPgPort, dbname, kPgPassword);
//     PGconn* conn = PQconnectdb(conninfo.c_str());
//     assert(PQstatus(conn) == CONNECTION_OK);
//     PGresult *res = NULL;
//     int seg = omp_get_thread_num();
//     long long start_count = seg * chunk + 1;
//     long long end_count = min((seg + 1) * chunk, N);
//     long long in_memory_count = ceilDiv(end_count-start_count+1, kInMemorySize);
//     for (long long i = 0; i < in_memory_count; i ++){
//       long long left = start_count + i * kInMemorySize;
//       long long right = min(left + kInMemorySize - 1, end_count);
//       sql = fmt::format("INSERT INTO {} SELECT * FROM {}({}, {});", Synthetic::table_name, kGenerateFunc, left, right);
//       res = PQexec(conn, sql.c_str());
//       PQclear(res);
//     }
//     PQfinish(conn);
//   }
//   pro.stop(2);
//   pro.clock(3);
//   {
//     string sql;
//     string conninfo = fmt::format("postgresql://{}@{}?port={}&dbname={}&password={}", kPgUser, kPgHostaddr, kPgPort, dbname, kPgPassword);
//     PGconn* conn = PQconnectdb(conninfo.c_str());
//     assert(PQstatus(conn) == CONNECTION_OK);
//     PGresult *res = NULL;

//     sql = fmt::format("ALTER TABLE {} ADD PRIMARY KEY ({});", Synthetic::table_name, kId);
//     res = PQexec(conn, sql.c_str());
//     PQclear(res);

//     PQfinish(conn);
//   }
//   pro.stop(3);
// }

void Synthetic::init(){
  // string sql;
  // string conninfo = fmt::format("postgresql://{}@{}?port={}&dbname={}&password={}", kPgUser, kPgHostaddr, kPgPort, dbname, kPgPassword);
  // PGconn* conn = PQconnectdb(conninfo.c_str());
  // assert(PQstatus(conn) == CONNECTION_OK);
  // PGresult *res = NULL;
  
  // sql = fmt::format("SET client_min_messages = warning;");
  // res = PQexec(conn, sql.c_str());
  // PQclear(res);

  // sql = "CREATE EXTENSION IF NOT EXISTS tablefunc;";
  // res = PQexec(conn, sql.c_str());
  // PQclear(res);

  // PQfinish(conn);
}