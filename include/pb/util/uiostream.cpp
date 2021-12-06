#define FMT_HEADER_ONLY

#include <filesystem>
#include <fstream>

#include "Eigen/Dense"
#include "uiostream.h"
#include "uconfig.h"
#include "fmt/core.h"
#include "ilcplex/ilocplex.h"
#include "umisc.h"
#include "fitsio.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;

/**
 * \brief   Return the filenames of all files that have the specified extension
 *          in the specified directory and all subdirectories.
 */
vector<string> getAllFiles(string root, string ext){
  vector<string> paths;
  for (auto &p : std::filesystem::recursive_directory_iterator(root)){
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
    problems[path] = (int) sol.size();
  }
  return problems;
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
  std::ifstream infile(fmt::format("{}/resource/solutions/{}.sol", kProjectHome, problem));
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

int getNumRows(string file_name){
  fitsfile *fptr;
  int status = 0;
  int hdunum, hdutype;
  long nrows;
  assert(!fits_open_file(&fptr, file_name.c_str(), READONLY, &status));
  if (fits_get_hdu_num(fptr, &hdunum) == 1) fits_movabs_hdu(fptr, 2, &hdutype, &status);
  else fits_get_hdu_type(fptr, &hdutype, &status);
  assert(hdutype != IMAGE_HDU);
  fits_get_num_rows(fptr, &nrows, &status);
  assert(!fits_close_file(fptr, &status));
  return nrows;
}

MatrixXd readTable(string file_name, const vector<int> &cols, vector<string> &column_names){
  fitsfile *fptr;
  int status = 0;
  int hdunum, hdutype, ncols, anynul;
  long nrows;
  assert(!fits_open_file(&fptr, file_name.c_str(), READONLY, &status));
  if (fits_get_hdu_num(fptr, &hdunum) == 1) fits_movabs_hdu(fptr, 2, &hdutype, &status);
  else fits_get_hdu_type(fptr, &hdutype, &status);
  assert(hdutype != IMAGE_HDU);
  fits_get_num_rows(fptr, &nrows, &status);
  fits_get_num_cols(fptr, &ncols, &status);
  MatrixXd table (nrows, cols.size()); // Column major
  char colname[100];
  char tpl[] = "*";
  int colnum = -1;
  int cur = 0;
  column_names.clear();
  while (status != COL_NOT_FOUND && cur < (int)cols.size()){
    fits_get_colname(fptr, CASEINSEN, tpl, colname, &colnum, &status);
    if (colnum == cols[cur]){
      string name (colname);
      column_names.push_back(name);
      cur ++;
    }
  }
  status = 0;
  for (int i = 0; i < (int)cols.size(); i ++){
    assert(!ffgcvd(fptr, cols[i]+1, 1, 1, nrows, -DBL_MAX, &table(0, i), &anynul, &status));
  }
  assert(!fits_close_file(fptr, &status));
  return table;
}