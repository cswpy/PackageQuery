#include "det_prob.h"
#include "pb/util/uconfig.h"
#include "pb/util/upostgres.h"
#include "pb/util/unumeric.h"

DetProb::~DetProb(){
}

DetProb::DetProb(){
}

DetProb::DetProb(int m, int n){
  resize(m, n);
}

void DetProb::resize(int m, int n){
  A.resize(m, n);
  bl.resize(m); bl.fill(-DBL_MAX);
  bu.resize(m); bu.fill(DBL_MAX);
  c.resize(n); c.fill(0);
  l.resize(n); l.fill(0);
  u.resize(n); u.fill(1);
  ids.resize(n);
}

void DetProb::uniformGenerate(int n, int expected_n, double att_var, double outlier_prob, bool restrict_count, bool is_positive, bool is_translate, int seed){
  if (restrict_count) resize(4, n);
  else resize(3, n);
  double left = -sqrt(3*att_var);
  double right = sqrt(3*att_var);
  if (is_positive){
    double translation = 0;
    if (is_translate) translation = 1;
    left = translation * right;
    right = (2 + translation) * right;
  }
  unsigned int local_seed;
  if (seed < 0){
    random_device rd;
    local_seed = rd();
  } else{
    local_seed = seed;
  }
  default_random_engine gen (local_seed);
  uniform_real_distribution dist(left, right);
  for (int i = 0; i < n; i ++){
    A(0, i) = dist(gen);
    A(1, i) = dist(gen);
    A(2, i) = dist(gen);
    if (restrict_count) A(3, i) = 1;
    c(i) = dist(gen);
  }
  normal_distribution bound_dist((left+right)/2*expected_n, sqrt(att_var*expected_n));
  //uniform_real_distribution dist(left*sqrt(expected_n), right*sqrt(expected_n));
  double tol = sqrt(att_var*expected_n/outlier_prob);
  double v = bound_dist(gen);
  bl(0) = v-tol;
  bu(0) = v+tol;
  bu(1) = bound_dist(gen)+tol;
  bl(2) = bound_dist(gen)-tol;
  if (restrict_count){
    bl(3) = expected_n/2;
    bu(3) = expected_n*3/2;
  }
}

void DetProb::normalGenerate(int n, int expected_n, double att_var, double outlier_prob, bool restrict_count, int seed){
  if (restrict_count) resize(4, n);
  else resize(3, n);
  unsigned int local_seed;
  if (seed < 0){
    random_device rd;
    local_seed = rd();
  } else{
    local_seed = seed;
  }
  default_random_engine gen (local_seed);
  normal_distribution dist(0.0, sqrt(att_var));
  for (int i = 0; i < n; i ++){
    A(0, i) = dist(gen);
    A(1, i) = dist(gen);
    A(2, i) = dist(gen);
    if (restrict_count) A(3, i) = 1;
    c(i) = dist(gen);
  }
  normal_distribution bound_dist(0.0, sqrt(expected_n*att_var));
  double tol = sqrt(att_var*expected_n/outlier_prob);
  double v = bound_dist(gen);
  bl(0) = v-tol;
  bu(0) = v+tol;
  bu(1) = bound_dist(gen)+tol;
  bl(2) = bound_dist(gen)-tol;
  if (restrict_count){
    bl(3) = expected_n/2;
    bu(3) = expected_n*3/2;
  }
}

// First column is the objective
void DetProb::tableGenerate(string table_name, vector<string>& cols, bool is_maximize, int n, int seed){
  PgManager pg = PgManager();
  long long size = pg.getSize(table_name);
  long long chunk = ceilDiv(size, (long long) kPCore);
  double probability = n / (double) size;
  unsigned int local_seed;
  if (seed < 0){
    random_device rd;
    local_seed = rd();
  } else{
    local_seed = seed;
  }
  seed_seq seq{local_seed};
  vector<unsigned int> local_seeds (kPCore);
  seq.generate(local_seeds.begin(), local_seeds.end());
  string col_name = join(cols, ",");
  int global_size = 0;
  MeanVar mv = MeanVar(3);
  #pragma omp parallel num_threads(kPCore)
  {
    int seg = omp_get_thread_num();
    default_random_engine gen (local_seeds[seg]);
    uniform_real_distribution dist (0.0, 1.0);
    long long start_id = seg * chunk + 1;
    long long end_id = (seg + 1) * chunk;
    vector<long long> selected_ids; selected_ids.reserve((int)(chunk * probability));
    for (long long id = start_id; id <= end_id; id ++){
      if (dist(gen) <= probability) selected_ids.push_back(id);
    }
    int start_index;
    #pragma omp critical (c1)
    {
      start_index = global_size;
      global_size += selected_ids.size();
    }
    #pragma omp barrier
    #pragma omp master
    {
      resize(4, global_size);
    }
    #pragma omp barrier
    vector<string> sids;
    for (auto id : selected_ids) sids.push_back(to_string(id));
    string sql = fmt::format("SELECT {},{} FROM \"{}\" WHERE {} IN ({})", kId, col_name, table_name, kId, join(sids, ","));
    PGconn *conn = PQconnectdb(pg.conninfo.c_str());
    assert(PQstatus(conn) == CONNECTION_OK);
    PGresult *res = NULL;
    res = PQexec(conn, sql.c_str());

    MeanVar local_mv = MeanVar(3);
    for (int i = 0; i < PQntuples(res); i++){
      int index = i + start_index;
      ids[index] = atol(PQgetvalue(res, i, 0));
      if (is_maximize) c(index) = atof(PQgetvalue(res, i, 1));
      else c(index) = -atof(PQgetvalue(res, i, 1));
      VectorXd x (3);
      for (int j = 0; j < 3; j ++){
        A(j, index) = atof(PQgetvalue(res, i, j+2));
        x(j) = A(j, index);
      }
      local_mv.add(x);
      A(3, index) = 1;
    }
    PQclear(res);
    PQfinish(conn);
    #pragma omp critical (c2)
    {
      mv.add(local_mv);
    }
  }
  vector<int> consSense = {LowerBounded, UpperBounded, Bounded};
  detBound = DetBound(consSense, mv.getMean(), mv.getVar(), seed);
}

double DetProb::boundGenerate(double E, double alpha, double hardness){
  bl(3) = E/2;
  bu(3) = E*3/2;
  return detBound.sampleHardness(E, alpha, hardness, bl, bu);
}

void DetProb::normalizeObjective(){
  double max_abs = 0;
  for (int i = 0; i < c.size(); i ++) max_abs = max(fabs(c(i)), max_abs);
  if (max_abs > 0) c /= max_abs;
}

void DetProb::truncate(){
  int n = (int) c.size();
  if (!n) return;
  int m = (int) bl.size();
  int new_m = 0;
  for (int i = 0; i < m; i ++){
    if (bl(i) != -DBL_MAX || bu(i) != DBL_MAX) new_m ++;
  }
  if (new_m < m){
    RMatrixXd new_A; new_A.resize(new_m, n);
    VectorXd new_bl (new_m);
    VectorXd new_bu (new_m);
    int index = 0;
    for (int i = 0; i < m; i ++){
      if (bl(i) != -DBL_MAX || bu(i) != DBL_MAX){
        memcpy(&(new_A(index, 0)), &(A(i, 0)), n*sizeof(double));
        new_bl(index) = bl(i);
        new_bu(index) = bu(i);
        index ++;
      }
    }
    A = new_A;
    bl = new_bl;
    bu = new_bu;
  } 
}

void DetProb::setSeed(int seed){
  detBound.setSeed(seed);
}
