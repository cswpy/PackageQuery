#include "pb/util/udeclare.h"
#include "pb/util/uconfig.h"
#include "pb/util/udebug.h"
#include "pb/det/kd_tree.h"
#include "pb/core/dual_reducer.h"
#include "pb/core/gurobi_solver.h"

using namespace pb;

int main(){
    KDTree kt;
    vector<string> cols = {"tmass_prox", "j", "h", "k"};
    //kt.partitionTable("ssds_demo", "s10_d5", cols, 10, 5);
    kt.partitionTable("ssds_6_1", "s5000_d10", cols, 5000, 10); 
    // Query: 0.1s, Build: 0.4s, I/O: 10s After refactoring: Query: 1.1s, Build: 0.4s, I/O: 0.6s
    //kt.partitionTable("ssds_7_1", "s20000_d20", cols, 20000, 20); 
    // Query: 10s, Build: 0.9s, I/O: 125s After refactoring: Query: 10.9s, Build: 7.9s, I/O: 5.7s
}