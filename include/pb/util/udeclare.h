#define FMT_HEADER_ONLY

#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <utility>
#include <string>
#include <numeric>
#include <algorithm>
#include <map>
#include <random>
#include <float.h>
#include <queue>
#include <omp.h>
#include <chrono>
#include <thread>
#include <ctime>
#include <cfloat>
#include <stdexcept>

#include "Eigen/Dense"
#include "fmt/core.h"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RMatrixXd;

enum ConsSense { LowerBounded, UpperBounded, Bounded};
enum SolStatus { NotFound, Found, Feasible, Infeasible, Unbounded, DualUnbounded, Timeout, NoPartitionFound, IncompatiblePartition};
enum FeasStatus { Unsolved, Feasibility, Infeasibility, LbConstraint, UbConstraint, LbVariable, UbVariable, Integrality, BadFilter};
enum GroupStatus { Unitialized, Unlocked, Locked};

namespace pb{
  using std::map;
  using std::unordered_map;
  using std::unordered_set;
  using std::priority_queue;
  using std::vector;
  using std::string;
  using std::pair;
  using std::tuple;
  using std::tie;
  using std::shared_ptr;
  
  using std::to_string;
  using std::iota;
  using std::sort;
  using std::min;
  using std::max;
  using std::map;
  using std::nth_element;
  using std::distance;
  using std::make_shared;

  using Eigen::VectorXd;
  using Eigen::VectorXi;
  using Eigen::MatrixXd;

  using std::random_device;
  using std::default_random_engine;
  using std::uniform_real_distribution;
  using std::normal_distribution;
  using std::seed_seq;

  using std::fstream;
  using std::ios;
  using std::cout;
  using std::cin;
  using std::endl;

  using std::invalid_argument;
}