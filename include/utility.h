#pragma once

#include <iostream>
#include <chrono>
#include <Eigen/Dense>
#include <thread>

using namespace Eigen;
using namespace std;

#define CORE_COUNT 80

constexpr char CPX_EQUAL = 'E';
constexpr char CPX_LESS_EQUAL = 'L';
constexpr char CPX_GREATER_EQUAL = 'G';
constexpr char CPX_RANGE = 'R';
constexpr double kFloatEps = 1e-10;
constexpr double kFeasTolerance = 1e-6;
constexpr double kStubbornMultiplier = 2;
constexpr double kSizeTol = 1e5;
constexpr int kRetargetUB = 1 << 10;
constexpr int kLatticeUB = 1 << 6;
constexpr int kBoundTol = 1;
constexpr int LS_NOT_FOUND = 0;
constexpr int LS_FOUND = 1;
constexpr int LS_INFEASIBLE = 2;
constexpr int LS_UNBOUNDED = 3;
constexpr int LS_FEASIBLE = 4;
constexpr bool kIgnoreDegenerate = true;
constexpr bool kDoPresolve = false;
const string kProjectHome = "/home/alm818/package_query";

bool isEqual(double x, double y);
bool isLess(double x, double y);
bool isGreater(double x, double y);
bool isLessEqual(double x, double y);
bool isGreaterEqual(double x, double y);

VectorXd readSolution(string problem);
unordered_map<string, int> getProblemSizes();
VectorXd bucketSort(vector<double> array);
string solMessage(int sol_status);

template <typename T>
double sign(T value) {
  if (value < 0) return -1.0;
  if (value > 0) return 1.0;
  return 0.0;
}

template <typename T>
string print(Matrix<T, Dynamic, 1> v) {
  string s = "[";
  for (int i = 0; i < v.size() - 1; i++) {
    s += to_string(v(i)) + " ";
  }
  if (v.size() > 0) {
    s += to_string(v(v.size() - 1));
  }
  s += "]\n";
  cout << s;
  return s;
}

template<typename F, typename... Args>
double exeTime(F func, Args&&... args) {
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  func(forward<Args>(args)...);
  return chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - t1).count() / 1000000.0;
}

