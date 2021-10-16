#define FMT_HEADER_ONLY

#include <sstream>
#include <string>
#include <fstream>
#include <random>
#include <limits.h>
#include <unistd.h>
#include <filesystem>

#include "ilcplex/ilocplex.h"
#include "utility.h"
#include <fmt/core.h>

namespace fs = filesystem;

bool isEqual(double x, double y, double eps){
  return fabs(x-y) < eps;
}

bool isLessEqual(double x, double y, double eps){
  return x-y < eps;
}

bool isGreaterEqual(double x, double y, double eps){
  return x-y > -eps;
}

bool isLess(double x, double y, double eps){
  return x-y < -eps;
}

bool isGreater(double x, double y, double eps){
  return x-y > eps;
}

vector<string> sol_messages = {
  "NOT FOUND",
  "FOUND",
  "INFEASIBLE",
  "UNBOUNDED",
  "FEASIBLE"
};

string solMessage(int sol_status){
  return sol_messages[sol_status];
}

string solCombination(VectorXd sol){
  string combs = "";
  for (int i = 0; i < sol.size(); i ++){
    int count = sol(i);
    if (count == 1) combs += to_string(i) + " ";
    else if (count > 1) combs += fmt::format("{}x{} ", to_string(i), count);
  }
  return combs;
}

const char* space = " \t\n\r\f\v";

// trim from end of string (right)
inline string& rtrim(string& s, const char* t = space){
  s.erase(s.find_last_not_of(t) + 1);
  return s;
}

// trim from beginning of string (left)
inline string& ltrim(string& s, const char* t = space){
  s.erase(0, s.find_first_not_of(t));
  return s;
}

// trim from both ends of string (right then left)
inline string& trim(string& s, const char* t = space){
  return ltrim(rtrim(s, t), t);
}

VectorXd readSolution(string problem) {
  int status = 0;
  CPXENVptr env = CPXXopenCPLEX(&status); assert(!status);
  CPXLPptr model = CPXXcreateprob(env, &status, problem.c_str()); assert(!status);
  assert(!CPXXreadcopyprob(env, model, fmt::format("{}/resource/test_cases/{}.mps", kProjectHome, problem).c_str(), NULL));
  int n = CPXXgetnumcols(env, model);
  long long surplus;
  assert(CPXXgetcolname(env, model, NULL, NULL, 0, &surplus, 0, n-1) == CPXERR_NEGATIVE_SURPLUS);
  long long storespace = -surplus;
  char** name = new char*[n];
  char* namestore = new char[storespace];
  assert(!CPXXgetcolname(env, model, name, namestore, storespace, &surplus, 0, n-1));
  assert(!CPXXfreeprob(env, &model));
  assert(!CPXXcloseCPLEX(&env));

  unordered_map<string, int> name_map;
  for (int i = 0; i < n; i ++){
    string var_name (name[i]);
    name_map[trim(var_name)] = i;
  }
  delete[] namestore;

  string line;
  VectorXd s (n); s.fill(0);
  ifstream infile(fmt::format("{}/resource/solutions/{}.sol", kProjectHome, problem));
  while (getline(infile, line)) {
    if (line.length() > 0 && line.at(0) != '='){
      auto pos = line.find_first_of(space);
      if (pos != string::npos){
        string val_str = line.substr(pos, string::npos);
        string name_str = line.substr(0, pos);
        double val = stod(trim(val_str));
        s(name_map[trim(name_str)]) = val;
      }
    }
  }
  infile.close();
  return s;
}

/**
 * \brief   Return the filenames of all files that have the specified extension
 *          in the specified directory and all subdirectories.
 */
vector<string> getAllFiles(string root, string ext){
  vector<string> paths;
  for (auto &p : fs::recursive_directory_iterator(root)){
    if (p.path().extension() == ext)
      paths.push_back(p.path().stem().string());
  }
  return paths;
}

unordered_map<string, int> getProblemSizes(){
  unordered_map<string, int> problems;
  auto paths = getAllFiles(fmt::format("{}/resource/solutions", kProjectHome), ".sol");
  for (auto path : paths){
    VectorXd sol = readSolution(path);
    problems[path] = sol.size();
  }
  return problems;
}

// Assuming uniform distribution [0,1]
// Return sorted array with 1 appended at the end
// O(n) on average
VectorXd bucketSort(vector<double> array){
  int n = array.size() + 1;
  VectorXd res (n); res(n-1) = 1.0;
  vector<vector<double>> buckets (n);
  for (int i = 0; i < (int)array.size(); i ++){
    int bucket_index = min((int)floor(array[i]*n), n-1);
    buckets[bucket_index].push_back(array[i]);
  }
  int index = 0;
  for (int i = 0; i < n; i ++){
    for (auto v : buckets[i]){
      res(index) = v;
      index ++;
    }
  }
  // Insertion sort
  for (int i = 1; i < n - 1; i ++){
    double key = res(i);
    int j = i - 1;
    while (j >= 0 && res(j) > key){
      res(j+1) = res(j);
      j --;
    }
    res(j+1) = key;
  }
  return res;
}