#include "pb/util/udeclare.h"
#include "pb/det/synthetic.h"
#include "pb/det/dlv.h"

using namespace pb;

void create_table(){
  Synthetic syn = Synthetic();
  long long N = 1000000000;
  syn.createMixed(N, 2, 2, 20, 100);
  syn.pro.print();
} 

void test_dlv(){
  DynamicLowVariance dlv = DynamicLowVariance();
  dlv.partition(Synthetic::table_name, "P0");
  dlv.pro.print();
}

int main(){
  //create_table();
  test_dlv();
}