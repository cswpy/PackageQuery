#include "pb/util/udeclare.h"
#include "pb/util/uconfig.h"
#include "pb/util/udebug.h"
#include "pb/det/synthetic.h"
#include "pb/core/dual_reducer.h"
#include "pb/core/gurobi_solver.h"

using namespace pb;

void create_table() {
  Synthetic syn = Synthetic();
  // vector<string> cols = {"tmass_prox", "j", "h", "k"};
  // syn.createSubtable("ssds", 6, cols, 1);

  vector<string> cols = {"price", "quantity", "discount", "tax"};
  syn.createSubtable("tpch", 3, cols, 1);
  syn.pro.print();
}

int main(){
  create_table();
}