#define FMT_HEADER_ONLY

#include <math.h>
#include <omp.h>
#include <regex>
#include <cfloat>
#include <string>
#include "fmt/core.h"
#include "umisc.h"
#include "boost/algorithm/string.hpp"

using boost::split;
using boost::is_any_of;
using std::stod;
using std::regex;
using std::smatch;
using std::to_string;

const string kPgDelim = ",";

string infAlias(double val, int precision){
  if (val == -DBL_MAX) return "-infinity";
  if (val == DBL_MAX) return "infinity";
  if (precision < 0) return to_string(val);
  return fmt::format("{:.{}Lf}", val, precision);
}

string join(vector<string> names, string delim){
  string res = "";
  if (names.size() == 0) return res;
  for (int i = 0; i < (int)(names.size()-1); i++) res += names[i] + delim;
  res += names[names.size()-1];
  return res;
}

string join(VectorXd vals, int precision){
  vector<string> str_vals;
  for (double v : vals) str_vals.push_back(fmt::format("{:.{}Lf}", v, precision));
  return join(str_vals, kPgDelim);
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

string pgJoin(vector<string> names){
  return fmt::format("'{{{}}}'", join(names, kPgDelim));
}

string pgJoin(VectorXd vals, int precision){
  return fmt::format("'{{{}}}'", join(vals, precision));
}

vector<string> pgStringSplit(char *s){
  string sub (s+1, strlen(s)-2);
  vector<string> res;
  split(res, sub, is_any_of(kPgDelim));
  return res;
}

VectorXd pgValueSplit(char *s){
  string sub (s+1, strlen(s)-2);
  vector<string> res;
  split(res, sub, is_any_of(kPgDelim));
  VectorXd vals (res.size());
  for (int i = 0; i < (int) res.size(); i ++) vals(i) = stod(res[i]);
  return vals;
}

pair<double, double> atop(char *s){
  VectorXd vals = pgValueSplit(s);
  double left = -DBL_MAX;
  if (!isinf(vals(0))) left = vals(0);
  double right = DBL_MAX;
  if (!isinf(vals(1))) right = vals(1);
  return {left, right};
}

bool isPartitionName(string table_name){
  regex expr ("^\\[([0-9]+).+\\](.+)");
  smatch m;
  regex_search(table_name, m, expr);
  return m.size() > 0;
}

string nextName(string table_name, string symbol){
  regex expr ("^\\[([0-9]+).+\\](.+)");
  smatch m;
  regex_search(table_name, m, expr);
  if (m.size() == 0) return fmt::format("[1{}]", symbol) + table_name;
  else{
    int index = stoi(m[1]) + 1;
    string name (m[2]);
    return fmt::format("[{}{}]{}", index, symbol, name);
  }
}

string nextGName(string table_name){
  return nextName(table_name, "G");
}

string nextPName(string table_name){
  return nextName(table_name, "P");
}

int getLayerIndex(string table_name){
  regex expr ("^\\[([0-9]+).+\\].+");
  smatch m;
  regex_search(table_name, m, expr);
  if (m.size() == 0) return 0;
  else return stoi(m[1]);
}

bool startsWith(string full, string sub){
  if (full.length() >= sub.length()){
    return full.compare(0, sub.length(), sub) == 0;
  }
  return false;
}

bool endsWith(string full, string sub){
  if (full.length() >= sub.length()){
    return full.compare(full.length() - sub.length(), sub.length(), sub) == 0;
  }
  return false;
}

bool isIn(vector<string> arr, string s){
  for (string a : arr){
    if (a == s) return true;
  }
  return false;
}

void allocate(char **vals, int index, const char *s){
  vals[index] = new char[strlen(s)];
  strcpy(vals[index], s);
}

void free(char** vals, int length){
  for (int i = 0; i < length; i ++) delete vals[i];
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