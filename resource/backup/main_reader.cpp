#define FMT_HEADER_ONLY

#include <iostream>
#include <Eigen/Dense>
#include <iostream>
#include <float.h>
#include <queue>
#include <random>
#include <chrono>
#include <typeinfo>
#include <numeric>

#include "map_sort.h"
#include "utility.h"
#include "parallel_pq.h"
#include "fmt/core.h"
#include "fitsio.h"
#include "omp.h"
#include "pseudo_walker.h"
#include "gurobi_solver.h"
#include "dual.h"
#include "dual_reducer.h"
#include "config.h"
#include "libpq-fe.h"
#include "problem.h"

using namespace std;
using namespace Eigen;

MeanVar readDataFromPG(string dbname, string table_name, string obj_att, vector<string> constr_atts, RMatrixXd& A, VectorXd& c, int n, VectorXi& indices){
  int nn;
  {
    string conninfo = fmt::format("postgresql://{}@{}?port={}&dbname={}&password={}", user, hostaddr, port, dbname, password);
    PGconn* conn = PQconnectdb(conninfo.c_str());
    assert(PQstatus(conn) == CONNECTION_OK);
    PGresult *res = NULL;
    string count_sql = fmt::format("SELECT COUNT(*) FROM {}", table_name);
    res = PQexec(conn, count_sql.c_str());
    nn = atoi(PQgetvalue(res, 0, 0));
    PQclear(res);
    PQfinish(conn);
  }
  vector<int> index (nn); iota(index.begin(), index.end(), 1);
  unsigned seed = chrono::system_clock::now().time_since_epoch().count();
  shuffle(index.begin(), index.end(), default_random_engine(seed));
  indices.resize(n);
  int m = constr_atts.size();
  A.resize(m, n);
  c.resize(n);
  vector<bool> is_count (m, false);
  vector<string> cols;
  for (int i = 0; i < m; i ++){
    if (constr_atts[i].length() == 0) is_count[i] = true;
    else cols.push_back(constr_atts[i]);
  }
  cols.push_back(obj_att);
  string col_names = join(cols);
  int nS = 0;
  VectorXd mean (m); mean.fill(0);
  VectorXd M2 (m); M2.fill(0);
  #pragma omp parallel num_threads(db_physical_core)
  {
    int db_chunk = n / omp_get_num_threads();
    string conninfo = fmt::format("postgresql://{}@{}?port={}&dbname={}&password={}", user, hostaddr, port, dbname, password);
    PGconn* conn = PQconnectdb(conninfo.c_str());
    assert(PQstatus(conn) == CONNECTION_OK);
    MeanVar mv = MeanVar(m);
    string select_sql = fmt::format("SELECT {} FROM {} WHERE {}.id=$1", col_names, table_name, table_name);
    string name = "select";
    Oid oid_types[1] = {23};
    PGresult *res = PQprepare(conn, name.c_str(), select_sql.c_str(), 1, (const Oid*)oid_types);
    PGresult *exec = NULL;
    const char* param_values[1];
    int param_lens[1];
    int param_formats[1] = {0};
    #pragma omp for nowait
    for (int i = 0; i < n; i ++){
      indices(i) = index[i];
      param_values[0] = to_string(indices(i)).c_str();
      param_lens[0] = strlen(param_values[0]);
      exec = PQexecPrepared(conn, name.c_str(), 1, param_values, param_lens, param_formats, 0);
      int cnt = 0;
      cout << indices(i) << " ";
      // for (int j = 0; j < m; j ++){
      //   if (!is_count[j]){
      //     A(j, indices(i)-1) = atof(PQgetvalue(exec, 0, cnt));
      //     cnt ++;
      //   } else A(j, indices(i)-1) = 1;
      // }
      // c(indices(i)-1) = atof(PQgetvalue(exec, 0, cols.size()-1));
      // mv.add(A.col(indices(i)-1));
      PQclear(exec);
    }
    PQclear(res);
    PQfinish(conn);
    #pragma omp critical
    {
      int tmp_nS = nS + mv.sample_count;
      VectorXd delta = mv.mean - mean;
      mean = (nS*mean + mv.sample_count*mv.mean) / tmp_nS;
      M2 += mv.M2 + delta.cwiseProduct(delta) * (nS * mv.sample_count / (double) (tmp_nS));
      nS = tmp_nS;
    }
  }
  MeanVar mv = MeanVar();
  mv.sample_count = nS;
  mv.M2 = M2;
  mv.mean = mean;
  mv.var = mv.M2 / (nS - 1);
  return mv;
}

void generateQuery(string dbname, string table_name, vector<string> constr_atts, string obj_att, vector<int> bounds, int k, RMatrixXd& A, VectorXd& bl, VectorXd& bu, VectorXd& c, int n, VectorXi& indices){
  MeanVar mv = readDataFromPG(dbname, table_name, obj_att, constr_atts, A, c, n, indices);
  // int m = A.outerSize();
  // int n = A.innerSize();
  // bl.resize(m);
  // bu.resize(m);
  // VectorXd mean = mv.mean * k;
  // VectorXd var = mv.var * k;
  // double outlier_prob = 0.6;
  // VectorXd tol = var*sqrt(1/outlier_prob);
  // for (int i = 0; i < m; i ++){
  //   if (constr_atts[i].length() > 0){
  //     if (bounds[i] <= 0) bl(i) = mean(i) - tol(i);
  //     else bl(i) = -DBL_MAX;
  //     if (bounds[i] >= 0) bu(i) = mean(i) + tol(i);
  //     else bu(i) = DBL_MAX;
  //   } else{
  //     if (bounds[i] <= 0) bl(i) = k*1.0/2;
  //     else bl(i) = -DBL_MAX;
  //     if (bounds[i] >= 0) bu(i) = k*3.0/2;
  //     else bu(i) = DBL_MAX;
  //   }
  // }
}

void testReader(){
  string dbname = "synthetic";
  string table_name = "n10000000_t1_v100";
  vector<string> constr_atts = {"a2", "a3", "a4", ""};
  vector<int> bounds = {0, -1, 1, 0};
  string obj_att = "a1";
  RMatrixXd A;
  VectorXd c, bl, bu;
  VectorXi indices;
  int n = 1000000;
  generateQuery(dbname, table_name, constr_atts, obj_att, bounds, 20, A, bl, bu, c, n, indices);
  // int m = A.outerSize();
  // VectorXd l (n); l.fill(0);
  // VectorXd u (n); u.fill(1);
  // cout << m << " " << n << endl;
  // // cout << A << endl;
  // print(bl);
  // print(bu);
  // DualReducer dr = DualReducer(db_physical_core, &A, &bl, &bu, &c, &l, &u, 1);
  // double r0_obj = dr.duals[0]->score;
  // double gap = (r0_obj - dr.best_score) / r0_obj * 100;
  // cout << solMessage(dr.status) << endl;
  // cout << solCombination(dr.best_x) << endl;
  // cout << dr.best_score << " " << r0_obj << " " << gap << "%" << endl;
  // GurobiSolver gs = GurobiSolver(A, bl, bu, c, l, u);
  // gs.solveRelaxed();
  // cout << gs.exe_relaxed << " " << gs.iteration_count << endl;
  // cout << gs.relaxed_cscore << endl;
}

int main(){
  testReader();
}