#include "pb/util/udeclare.h"
#include "pb/det/synthetic.h"
#include "pb/det/dlv.h"

using namespace pb;

void create_table(){
  string dbname = "benchmark";
  Synthetic syn = Synthetic(dbname);
  long long N = 1000000;
  double var = 100;
  syn.createUniform(N, sqrt(3*var), var);
  syn.pro.print();
} 

void test_dlv(){
  string dbname = "benchmark";
  DynamicLowVariance dlv = DynamicLowVariance(dbname);
  dlv.pro.print();
}

int main(){
  // create_table();
  test_dlv();
}