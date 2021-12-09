#define FMT_HEADER_ONLY

#include "fmt/core.h"
#include "dlv.h"
#include "pb/util/uconfig.h"
#include "pb/util/umisc.h"
#include "pb/util/unumeric.h"
#include "pb/util/udeclare.h"

using namespace pb;

double DynamicLowVariance::kGroupRatio = 0.01;
int DynamicLowVariance::kMinPartitionSize = 1000000;
int DynamicLowVariance::kInitialSize = 100000;

DynamicLowVariance::~DynamicLowVariance(){
  PQfinish(_conn);
  delete pg;
}

DynamicLowVariance::DynamicLowVariance(string dbname): dbname(dbname){
  vector<string> names = {"Init", "Partition", "2", "3", "4", "5", "6"};
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

  _sql = fmt::format(""
    "CREATE TABLE IF NOT EXISTS {}("
    "	table_name VARCHAR(63) NOT NULL,"
    "	count BIGINT DEFAULT 0,"
    " column_names text[],"
    "	mean DOUBLE PRECISION[],"
    "	M2 DOUBLE PRECISION[]"
    ")", kStatTable);
  _res = PQexec(_conn, _sql.c_str());
  PQclear(_res);

  _sql = fmt::format(""
    "CREATE TABLE IF NOT EXISTS {}("
    "	table_name VARCHAR(63) NOT NULL,"
    "	partition_name VARCHAR(31) NOT NULL,"
    "	column_names text[],"
    "	group_ratio DOUBLE PRECISION,"
    "	min_partition_size INTEGER,"
    "	initial_size INTEGER,"
    " layer_count INTEGER"
    ");", kPartitionTable);
  _res = PQexec(_conn, _sql.c_str());
  PQclear(_res);
  pro.stop(0);
}

void DynamicLowVariance::writeStats(string table_name, Stat *stat){
  
}

void DynamicLowVariance::partition(string table_name){
  pro.clock(1);
  Stat *stat = pg->computeStats(table_name);
  print(stat->mean);
  print(stat->M2);
  writeStats(table_name, stat);
  pro.stop(1);
}