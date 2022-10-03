#include "det_sql.h"

#include "pb/util/uconfig.h"
#include "pb/util/umisc.h"

DetSql::~DetSql(){
  delete pg;
}

DetSql::DetSql(string table_name, string obj_col, bool is_maximize, vector<string> att_cols, vector<int> att_senses, bool has_count_constraint, long long u)
  : table_name(table_name), obj_col(obj_col), is_maximize(is_maximize), att_cols(att_cols), att_senses(att_senses), u(u), has_count_constraint(has_count_constraint){
  pg = new PgManager();
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
  
bool DetSql::isFiltering(){
  return filter_cols.size() > 0;
}