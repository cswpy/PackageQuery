#include "uconfig.h"

const double kPostgresMemoryScale = 8.0;

long long getTupleCount(int tuple_size, int core_count, double main_memory){
  return (long long) (main_memory / (double) (core_count * tuple_size * 8 * kPostgresMemoryScale) * 1000000000);
}