#pragma once

#include <vector>
#include <string>
#include <unordered_map>

#include "Eigen/Dense"

using std::vector;
using std::string;
using std::unordered_map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Cplex
vector<string> getAllFiles(string root, string ext);
unordered_map<string, int> getProblemSizes();

//ccfits
MatrixXd readTable(string file_name, const vector<int> &cols, vector<string> &column_names);
int getNumRows(string file_name);
VectorXd readSolution(string problem);