#include "uconfig.h"

long long getTupleCount(int tuple_size, int core_count){
  return (long long) (kMainMemorySize / (double) (core_count * tuple_size * 8) * 1000000000);
}