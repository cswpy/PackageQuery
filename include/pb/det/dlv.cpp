#define FMT_HEADER_ONLY

#include "fmt/core.h"
#include "dlv.h"
#include "libpq-fe.h"
#include "pb/util/uconfig.h"

double DynamicLowVariance::kGroupRatio = 0.01;
int DynamicLowVariance::kMinPartitionSize = 1000000;
int DynamicLowVariance::kInitialSize = 100000;

const string kStatTable = "dlv_stats";
const string kPartitionTable = "dlv_partitions";

void DynamicLowVariance::init(){
  pro.clock(0);
  string sql;
  string conninfo = fmt::format("postgresql://{}@{}?port={}&dbname={}&password={}", kPgUser, kPgHostaddr, kPgPort, dbname, kPgPassword);
  PGconn* conn = PQconnectdb(conninfo.c_str());
  assert(PQstatus(conn) == CONNECTION_OK);
  PGresult *res = NULL;

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

  pro.stop(0);
}

DynamicLowVariance::DynamicLowVariance(string dbname): dbname(dbname){
  vector<string> names = {"Init", "1", "2", "3", "4", "5", "6"};
  pro = Profiler(names);
  init();
}