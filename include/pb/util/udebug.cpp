#define FMT_HEADER_ONLY

#include <omp.h>
#include <cmath>

#include "udebug.h"
#include "pb/util/unumeric.h"
#include "fmt/core.h"

using std::min;

static vector<string> sol_messages = {
  "NOT FOUND",
  "FOUND",
  "FEASIBLE",
  "INFEASIBLE",
  "UNBOUNDED",
  "DUAL UNBOUNDED",
  "TIMEOUT",
  "NO PARTITION FOUND",
  "INCOMPATIBLE PARTITION",
  "NUMERICAL UNSTABILITY",
  "LP FOUND BUT NO ILP FOUND"
};

static vector<string> feas_messages = {
  "Unsolved",
  "Feasibility",
  "Infeasibility",
  "Lower bound constraint violation",
  "Upper bound constraint violation",
  "Lower bound variable violation",
  "Upper bound constraint violation",
  "Integrality violation",
  "Filter condition violation"
};

void showHistogram(vector<int> x, int bucket_count, double start, double end){
  if (x.size() == 0) return;
  VectorXd dx (x.size());
  for (int i = 0; i < (int)x.size(); i ++) dx(i) = x[i];
  showHistogram(dx, bucket_count, start, end);
}

void showHistogram(VectorXd x, int bucket_count, double start, double end){
  if (x.size() == 0) return;
  assert(bucket_count >= 1);
  double min_x = x.minCoeff();
  double max_x = x.maxCoeff();
  if (start <= min_x) start = min_x;
  if (end >= max_x) end = max_x;
  if (end <= start){
    start = min_x;
    end = max_x;
  }
  int outlier = 0;
  vector<int> buckets (bucket_count, 0);
  double dist = (end-start) / bucket_count;
  for (double x_val : x){
    if (start <= x_val && x_val <= end){
      int index = (int) min(floor((x_val - start) / dist), bucket_count-1.0);
      buckets[index] ++;
    } else outlier ++;
  }
  for (int i = 0; i <= bucket_count; i ++){
    double val = start + i * dist;
    fmt::print("{: 8.2Lf}", val);
  }
  fmt::print("\n");
  for (int i = 0; i < bucket_count; i ++){
    fmt::print("{: 8d}", buckets[i]);
  }
  fmt::print("{: 8d}\n", outlier);
}

string solMessage(int sol_status){
  return sol_messages[sol_status];
}

string feasMessage(int feas_status){
  return feas_messages[feas_status];
}

string solCombination(VectorXd sol){
  string combs = "";
  for (int i = 0; i < sol.size(); i ++){
    int count = (int) round(sol(i));
    if (count == 1) combs += to_string(i) + " ";
    else if (count > 1) combs += fmt::format("{}x{} ", to_string(i), count);
  }
  return combs;
}

string showClassification(double x) {
  switch(std::fpclassify(x)) {
    case FP_INFINITE:  return "Inf";
    case FP_NAN:       return "NaN";
    case FP_NORMAL:    return "normal";
    case FP_SUBNORMAL: return "subnormal";
    case FP_ZERO:      return "zero";
    default:           return "unknown";
  }
}

void Profiler::init(int n){
  time.resize(n); time.fill(0);
  count.resize(n); count.fill(0);
  eps.clear(); eps.resize(n);
}

Profiler::Profiler(){
}

Profiler::Profiler(int n): n(n){
  names.resize(n);
  for (int i = 0; i < n; i ++) names[i] = to_string(i);
  init(n);
}

Profiler::Profiler(vector<string> names): names(names){
  n = (int) names.size();
  init(n);
}

void Profiler::clock(int i, bool is_parallel){
  if (is_parallel){
    #pragma omp barrier
    #pragma omp master
    {
      if (0 <= i && i < n) eps[i] = std::chrono::high_resolution_clock::now();
    }
    #pragma omp barrier
  } else{
    if (0 <= i && i < n) eps[i] = std::chrono::high_resolution_clock::now();
  }
}

void Profiler::stop(int i, bool is_parallel){
  if (is_parallel){
    #pragma omp barrier
    #pragma omp master
    {
      if (0 <= i && i < n){
        time(i) += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - eps[i]).count() / 1000000.0;
        count(i) ++;
      }
    }
    #pragma omp barrier
  } else{
    if (0 <= i && i < n){
      time(i) += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - eps[i]).count() / 1000000.0;
      count(i) ++;
    }
  }
}

void Profiler::print(){
  for (int i = 0; i < n; i ++){
    if (count(i) > 0){
      fmt::print("{}-{} Total:{:.4Lf}ms Count:{} Average:{:.4Lf}ms\n", i, names[i], time(i), count(i), time(i)/count(i));
    }
  }
}

void Profiler::add(Profiler &profiler, int core){
  assert(n == profiler.n);
  for (int i = 0; i < n; i ++){
    time(i) += profiler.time(i) / core;
    count(i) += ceilDiv(profiler.count(i), core);
  }
}

#ifdef _WIN32
  #include "windows.h"
  #include "psapi.h"

  double currentRAM(){
    PROCESS_MEMORY_COUNTERS_EX pmc;
    SIZE_T ram = pmc.WorkingSetSize; // In bytes
    return (double) ram / (1 << 30);
  }
#else // Assuming Linux
  double currentRAM(){
    FILE* file = fopen("/proc/self/status", "r");
    double result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
      if (strncmp(line, "VmRSS:", 6) == 0){
        string s = line;
        s = s.substr(7, s.find('k') - 7);
        trim(s); // In Kbytes
        result = stod(s) / (1 << 20);
        break;
      }
    }
    fclose(file);
    return result;
  }
#endif