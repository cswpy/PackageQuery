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
#include <pqxx/pqxx>

#include "map_sort.h"
#include "utility.h"
#include "parallel_pq.h"
#include "fmt/core.h"
#include "fitsio.h"
#include "omp.h"
#include "lattice_solver.h"
#include "pseudo_walker.h"
#include "simplex.h"
#include "gurobi_lattice_solver.h"
#include "gurobi_solver.h"
#include "dual.h"
#include "dual_reducer.h"
#include "reducer.h"
#include "fmt/core.h"
#include "config.h"

using namespace std;
using namespace Eigen;
using namespace pqxx;

void generateQuery(int qi, int n, RMatrixXd& A, VectorXd& bl, VectorXd& bu, VectorXd& c){
  double outlier_prob = 0.6;
  vector<double> lbs = {-1, -1, -1, -1, 0, 0, 0, 0};
  vector<int> expected_ns = {20, 20, 500, 500, 20, 20, 500, 500};
  vector<double> att_vars = {1, 100, 1, 100, 1, 100, 1, 100};
  int expected_n = expected_ns[qi];
  double att_var = att_vars[qi];
  A.resize(4, n); bl.resize(4); bu.resize(4); c.resize(n); 
  default_random_engine gen {static_cast<long unsigned int>(time(0))};
  double left = lbs[qi];
  double right = 1.0;
  uniform_real_distribution u_dist(left, right);
  double multiplicity = att_var / (right-left) * 12;
  //normal_distribution att_dist(10.0, att_var);
  int expected_numvar = expected_n;
  double mean = 0.5*multiplicity*expected_numvar;
  double var = att_var*expected_numvar;
  normal_distribution n_dist(0.0, var);
  //normal_distribution n_dist_c(0.0, 1.0/12);
  uniform_real_distribution n_dist_c(left, right);
  #pragma omp parallel for num_threads(CORE_COUNT)
  for (int i = 0; i < n; i ++){
    A(0, i) = u_dist(gen)*multiplicity;
    A(1, i) = u_dist(gen)*multiplicity;
    A(2, i) = u_dist(gen)*multiplicity;
    A(3, i) = 1;
    c(i) = n_dist_c(gen)*multiplicity;
    //c(i) = A(0, i) + A(1, i) + A(2, i);
  }
  double tol = var*sqrt(1/outlier_prob);
  bl(0) = mean - tol;
  bu(0) = mean + tol;
  bl(1) = -DBL_MAX;
  bu(1) = mean + tol;
  bl(2) = mean- tol;
  bu(2) = DBL_MAX;
  bl(3) = expected_n*1.0/2;
  bu(3) = expected_n*3.0/2;
}

void mapSort(){
  vector<string> names = {"1", "2"};
  Profiler pro = Profiler(names);
  default_random_engine gen {static_cast<long unsigned int>(time(0))};
  uniform_real_distribution u_dist(0.0, 1.0);
  int N = 100000000;
  int repeat = 1;
  map_sort::MapSort<pair<double,int>> map_sort;
  map_sort.Init(N);
  for (int j = 0; j < repeat; j ++){
    pair<double,int>* x = new pair<double,int>[N];
    for (int i = 0; i < N ; i ++) x[i] = {u_dist(gen), i};
    pair<double,int>* y = new pair<double,int>[N];
    for (int i = 0; i < N ; i ++) y[i] = {x[i].first, i};
    pro.clock(0, false);
    map_sort.Sort(x, N, 8);
    pro.stop(0, false);
    pro.clock(1, false);
    sort(y, y+N);
    pro.stop(1, false);
    for (int i = 0; i < N; i ++) assert(x[i].first == y[i].first);
    delete[] x;
    delete[] y;
  }
  pro.print();
}

void partitionTest(){
  default_random_engine gen {static_cast<long unsigned int>(time(0))};
  double var = 100;
  normal_distribution n_dist(0.0, var);
  int N = 100000;
  VectorXd a (N);
  for (int i = 0; i < N; i ++){
    a(i) = n_dist(gen);
  }
  double reduce_factor = 2;
  double lim_var = var/reduce_factor;
  vector<int> delim;
  sort(a.begin(), a.end());
  ScalarMeanVar mv = ScalarMeanVar();
  for (int i = 0; i < N; i ++){
    mv.add(a(i));
    if (mv.var > lim_var){
      mv.reset();
      delim.push_back(i);
    }
  }
  for (int i = 0; i < delim.size(); i ++) cout << delim[i] << " ";
  cout << endl;
  int cur = 0;
  mv.reset();
  for (int i = 0; i < delim.size(); i ++){
    while (cur < delim[i]){
      mv.add(a(cur));
      cur ++;
    }
    cout << i << " "<< mv.mean << " " << mv.var << endl;
    mv.reset();
  }
  showHistogram(a, 20, 0, 0);
}

void createSyntheticData(connection& C, int n, int type, double var){
  try{
    work W(C);
    string table_name = fmt::format("N{}_T{}_V{}", n, type, (int)var);

    string clear_sql = fmt::format("DROP TABLE IF EXISTS {};", table_name);
    W.exec(clear_sql);

    string create_sql = fmt::format("CREATE TABLE IF NOT EXISTS {} ("\
    "id SERIAL PRIMARY KEY,"\
    "a1 DOUBLE PRECISION,"\
    "a2 DOUBLE PRECISION,"\
    "a3 DOUBLE PRECISION,"\
    "a4 DOUBLE PRECISION"\
    ");", table_name);
    W.exec(create_sql);

    double left, right;
    if (type == 1){
      left = -sqrt(3*var);
      right = sqrt(3*var);
    } else if (type == 2){
      left = 0;
      right = 2*sqrt(3*var);
    }
    W.commit();

    default_random_engine gen {static_cast<long unsigned int>(time(0))};
    uniform_real_distribution u_dist(left, right);
    string insert_sql = fmt::format("INSERT INTO {}(a1,a2,a3,a4) VALUES ", table_name);
    int bulk = 1000000;
    for (int i = 0; i < n/bulk; i ++){
      work Wi (C);
      string bulk_sql = insert_sql;
      for (int j = 0; j < bulk-1; j ++){
        bulk_sql += fmt::format("({:.16Lf},{:.16Lf},{:.16Lf},{:.16Lf}),", u_dist(gen), u_dist(gen), u_dist(gen), u_dist(gen));
      }
      bulk_sql += fmt::format("({:.16Lf},{:.16Lf},{:.16Lf},{:.16Lf});", u_dist(gen), u_dist(gen), u_dist(gen), u_dist(gen));
      Wi.exec(bulk_sql);
      Wi.commit();
      fmt::print("{}%\n", 100.0*(i+1)*bulk/n);
    }
  } catch (const exception& e){
    cerr << e.what() << endl;
  }
}

void postgresql(){
  try{
    string dbname = "synthetic";
    connection C (fmt::format("dbname={} user={} password={} hostaddr={} port={}", dbname, user, password, hostaddr, port));
    assert(C.is_open());
    createSyntheticData(C, 100000000, 1, 1.0);
    createSyntheticData(C, 100000000, 1, 100.0);
    createSyntheticData(C, 100000000, 2, 1.0);
    createSyntheticData(C, 100000000, 2, 100.0);
  } catch (const exception& e){
    cerr << e.what() << endl;
  }
}

double kGroupTolerance = 10;

void dynamic_kd_partition(string dbname, string table_name, int n, vector<string> cols, double group_ratio, double var_ratio){
  vector<string> names = {"0", "1", "2", "3", "4", "5"};
  Profiler pro = Profiler(names);
  pro.clock(0);
  int chunk = n / db_logical_core;
  int m = cols.size();
  RMatrixXd A (m, n);
  int nS = 0;
  VectorXd mean (m); mean.fill(0);
  VectorXd M2 (m); M2.fill(0);
  #pragma omp parallel num_threads(db_logical_core)
  {
    try{
      connection C (fmt::format("dbname={} user={} password={} hostaddr={} port={}", dbname, user, password, hostaddr, port));
      assert(C.is_open());
      string col_names = "";
      for (int i = 0; i < cols.size()-1; i++) col_names += cols[i] + ",";
      col_names += cols[cols.size()-1];
      nontransaction N (C);
      MeanVar mv = MeanVar(m);
      #pragma omp for nowait
      for (int i = 0; i < db_logical_core; i ++){
        int left = i * chunk;
        int right = (i+1)*chunk;
        string select_sql = fmt::format("SELECT {} FROM {} WHERE id BETWEEN {} and {};", col_names, table_name, left+1, right);
        result R (N.exec(select_sql));
        for (int j = 0; j < R.size(); j ++){
          const auto& tup = R[j];
          for (int k = 0; k < m; k ++){
            string str_val (tup[k].view());
            A(k, j+left) = stod(str_val);
          }
          mv.add(A.col(j+left));
        }
      }
      #pragma omp critical
      {
        int tmp_nS = nS + mv.sample_count;
        VectorXd delta = mv.mean - mean;
        mean = (nS*mean + mv.sample_count*mv.mean) / tmp_nS;
        M2 += mv.M2 + delta.cwiseProduct(delta) * (nS * mv.sample_count / (double) (tmp_nS));
        nS = tmp_nS;
      }
    } catch (const exception& e){
      cerr << e.what() << endl;
    }
  }
  pro.stop(0);
  pro.clock(1);
  // Phase-1
  int partition_index = -1;
  double max_var = -1;
  for (int i = 0; i < m; i ++){
    double var = M2(i) / (nS - 1);
    if (max_var < var){
      max_var = var;
      partition_index = i;
    }
  }
  pair<double, int>* pis = new pair<double, int>[n];
  #pragma omp parallel for num_threads(db_physical_core)
  for (int i = 0; i < n; i ++) pis[i] = {A(partition_index, i), i};
  map_sort::Sort(pis, n, db_physical_core);
  pro.stop(1);
  int group_lim = (int) ceil((n * group_ratio) + db_physical_core * kGroupTolerance / var_ratio);
  int group_count = 0;
  vector<vector<double>> low (m, vector<double>(group_lim, -DBL_MAX));
  vector<vector<double>> high (m, vector<double>(group_lim, DBL_MAX));

  delete[] pis;
  pro.print();
  // Phase-2
}

struct IndexComp{
  const RMatrixXd& mat;
  int j;
  IndexComp(const RMatrixXd& mat, int j): mat(mat), j(j){
  }
  inline bool operator()(int i1, int i2){
    return mat(j, i1) < mat(j, i2);
  }
};

void customSort(){
  // VectorXd a (5); a << 3, 1, 0, 5, -1;
  // VectorXi b (5); b << 0,1,2,3,4;
  // print(b);
  // b.conservativeResize(3);
  // print(b);
  // sort(b.begin(), b.end(), comp(a));
  // print(b);
  // tuple<double, double, double> a = {3,4,3};
  // tuple<double, double, double> b = {3,2,3};
  // double c1,c2,c3;
  // tie(c1,c2,c3) = a;
  // cout << c1 << " " << c2 << " " << c3 << endl;
}

int main(){
  //partitionTest();
  postgresql();
  //mapSort();
  //vector<string> cols = {"a1", "a2", "a3", "a4"};
  //dynamic_kd_partition("synthetic", "n1000000_t1_v1", 1000000, cols, 0.01, 0.2);
  //customSort();
  //fmt::print("{:.16Lf}\n", -DBL_MAX);
}