#pragma once

#include <vector>
#include <string>
#include <chrono>
#include <utility>

#include "Eigen/Dense"

inline char separator(){
  #ifdef _WIN32
      return '\\';
  #else
      return '/';
  #endif
}

using std::vector;
using std::string;
using Eigen::VectorXd;

static const char* space = " \t\n\r\f\v";

string infAlias(double val, int precision);
string join(vector<string> names, string delim);
string join(VectorXd vals, int precision);
string& trim(string& s, const char *t=space);
string pgJoin(vector<string> names);
string pgJoin(VectorXd vals, int precision);
vector<string> pgStringSplit(char *s);
VectorXd pgValueSplit(char *s);
string nextGName(string table_name);
string nextPName(string table_name);
int getLayerIndex(string table_name);

bool startsWith(string full, string sub);
bool endsWith(string full, string sub);

// Return the number of threads that would be executed in parallel regions
int GetMaxThreads();
// Set the number of threads that would be executed in parallel regions
void SetNumThreads(int num_threads);
// Return the thread number, which lies in [0, the number of threads)
int GetThreadId();

template<typename F, typename... Args>
double exeTime(F func, Args&&... args) {
  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  func(std::forward<Args>(args)...);
  return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - t1).count() / 1000000.0;
}