#define FMT_HEADER_ONLY

#pragma once

#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <string>
#include <numeric>
#include <algorithm>
#include <random>
#include <float.h>
#include <omp.h>
#include <chrono>
#include <ctime>
#include <cfloat>

#include "Eigen/Dense"
#include "fmt/core.h"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RMatrixXd;

enum SolStatus { NotFound, Found, Feasible, Infeasible, Unbounded, DualUnbounded, Timeout};
enum FeasStatus { Unsolved, Feasibility, Infeasibility, LbConstraint, UbConstraint, LbVariable, UbVariable, Integrality};
enum GroupStatus { Unitialized, Unlocked, Locked};

namespace pb{
  using std::cout;
  using std::cin;
  using std::endl;

  using std::unordered_map;
  using std::unordered_set;
  using std::vector;
  using std::string;
  using std::pair;
  using std::tuple;
  using std::tie;
  
  using std::to_string;
  using std::iota;
  using std::sort;
  using std::min;
  using std::max;

  using Eigen::VectorXd;
  using Eigen::VectorXi;
  using Eigen::MatrixXd;
}