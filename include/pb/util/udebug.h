#pragma once

#define UNUSED(expr) (void)(expr)

#include <vector>
#include <string>
#include <chrono>
#include <iostream>
#include <cfloat>
#include <unordered_map>
#include <map>

#include "pb/util/umisc.h"
#include "pb/util/uconfig.h"

#include "Eigen/Dense"

#define NUM_CLOCK 10

#if DEBUG
  #define INIT_CLOCK(pro) pro = Profiler(NUM_CLOCK)
  #define CREATE_CLOCK(pro) Profiler pro = Profiler(NUM_CLOCK)
  #define START_CLOCK(pro, i) pro.clock(i)
  #define END_CLOCK(pro, i) pro.stop(i)
  #define ADD_CLOCK(local_pro, core)\
  ({                                \
    _Pragma("omp critical")         \
    {                               \
      pro.add(local_pro, core);     \
    }                               \
  })

#else
  #define INIT_CLOCK(pro) (void)0
  #define CREATE_CLOCK(pro) (void)0
  #define START_CLOCK(pro, i) (void)0
  #define END_CLOCK(pro, i) (void)0            
  #define ADD_CLOCK(local_pro, core) (void)0
#endif

using std::cout;
using std::vector;
using std::unordered_map;
using std::map;
using std::string;
using std::to_string;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::Matrix;
using Eigen::Dynamic;

void showHistogram(vector<int> x, int bucket_count, double start, double end);
void showHistogram(VectorXd x, int bucket_count, double start, double end);

string solMessage(int sol_status);
string feasMessage(int feas_status);
string solCombination(VectorXd sol);
string showClassification(double x);

template <typename T>
string print(vector<T> v) {
  string s = "[";
  for (int i = 0; i < (int) (v.size() - 1); i++)  s += infAlias((double) v[i], -1) + " ";
  if (v.size() > 0) s += infAlias((double) v[v.size() - 1], -1);
  s += "]\n";
  cout << s;
  return s;
}

template <typename T>
string print(Matrix<T, Dynamic, 1> v) {
  string s = "[";
  for (int i = 0; i < (int) (v.size() - 1); i++)  s += infAlias((double) v(i), -1) + " ";
  if (v.size() > 0) s += infAlias((double) v(v.size() - 1), -1);
  s += "]\n";
  cout << s;
  return s;
}

template <typename T>
string print(unordered_map<long long, T> sol){
  string s = "{ ";
  for (const auto& [key, value] : sol){
    s += to_string(key) + ":" + to_string(value) + " ";
  }
  s += "}\n";
  cout << s;
  return s;
}

template <typename T>
string print(map<long long, T> sol){
  string s = "{ ";
  for (const auto& [key, value] : sol){
    s += to_string(key) + ":" + to_string(value) + " ";
  }
  s += "}\n";
  cout << s;
  return s;
}

template <typename T>
string shortPrint(Matrix<T, Dynamic, 1> v) {
  string s = "[";
  for (int i = 0; i < v.size() - 1; i++) {
    if (isEqual(v(i), 0)) continue;
    s += to_string(i) + ":" + to_string(v(i)) + " ";
  }
  if (v.size() > 0 && !isEqual(v(v.size() - 1), 0)) {
    s += to_string(v.size() - 1) + ":" + to_string(v(v.size() - 1));
  }
  s += "]\n";
  cout << s;
  return s;
}

class Profiler{
public:
  int n;
  VectorXd time;
  VectorXi count;
  vector<std::chrono::high_resolution_clock::time_point> eps;
  vector<string> names;
private:
  void init(int n);
public:
  Profiler();
  Profiler(int n);
  Profiler(vector<string> names);
  void clock(int i, bool is_parallel=false);
  void stop(int i, bool is_parallel=false);
  void print();
  void add(Profiler &profiler, int core);
};

double currentRAM();
