#define FMT_HEADER_ONLY

#include <omp.h>
#include "fmt/core.h"
#include "umisc.h"

string join(vector<string> names, string delim){
  string res = "";
  for (int i = 0; i < names.size()-1; i++) res += names[i] + delim;
  res += names[names.size()-1];
  return res;
}

// trim from end of string (right)
string& rtrim(string &s, const char *t){
  s.erase(s.find_last_not_of(t) + 1);
  return s;
}

// trim from beginning of string (left)
string& ltrim(string& s, const char *t){
  s.erase(0, s.find_first_not_of(t));
  return s;
}

// trim from both ends of string (right then left)
string& trim(string& s, const char *t){
  return ltrim(rtrim(s, t), t);
}

string pgJoin(vector<string> names, string delim){
  return fmt::format("'{{{}}}'", join(names, delim));
}

string pgJoin(VectorXd vals, int precision){
  vector<string> str_vals;
  for (double v : vals) str_vals.push_back(fmt::format("{:.{}Lf}", v, precision));
  return fmt::format("'{{{}}}'", join(str_vals, ","));
}


// Return the number of threads that would be executed in parallel regions
int GetMaxThreads() {
  #ifdef _OPENMP
    return omp_get_max_threads();
  #else
    return 1;
  #endif
}

// Set the number of threads that would be executed in parallel regions
void SetNumThreads(int num_threads) {
  #ifdef _OPENMP
    omp_set_num_threads(num_threads);
  #else
    if (num_threads != 1) {
      assert(!"compile with -fopenmp");
    }
  #endif
}

// Return the thread number, which lies in [0, the number of threads)
int GetThreadId() {
  #ifdef _OPENMP
    return omp_get_thread_num();
  #else
    return 0;
  #endif
}