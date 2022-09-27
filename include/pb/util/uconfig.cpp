#include "uconfig.h"

const double kPostgresMemoryScale = 16.0;

long long getTupleCount(int tuple_size, int core_count){
  return (long long) (kMainMemorySize / (double) (core_count * tuple_size * 8 * kPostgresMemoryScale) * 1000000000);
}