#include "det_prob.h"

#include "pb/util/uconfig.h"
#include "pb/util/upostgres.h"
#include "pb/util/unumeric.h"

DetProb::~DetProb(){
}

DetProb::DetProb(){
}

void DetProb::operator=(const DetProb& dp){
  this->det_sql = dp.det_sql;
  this->A = dp.A;
  this->bl = dp.bl;
  this->bu = dp.bu;
  this->c = dp.c;
  this->l = dp.l;
  this->u = dp.u;
  this->ids = dp.ids;
  this->det_bound = dp.det_bound;
}

DetProb::DetProb(const DetProb& dp){
  this->det_sql = dp.det_sql;
  this->A = dp.A;
  this->bl = dp.bl;
  this->bu = dp.bu;
  this->c = dp.c;
  this->l = dp.l;
  this->u = dp.u;
  this->ids = dp.ids;
  this->det_bound = dp.det_bound;
}

DetProb::DetProb(int m, int n){
  resize(m, n);
}

DetProb::DetProb(DetSql &det_sql, long long n, int seed): det_sql(&det_sql){
  setSeed(seed);
  PgManager pg = PgManager();
  long long size = pg.getSize(det_sql.table_name);
  long long chunk = ceilDiv(size, (long long) kPCore);
  double probability = 1.0;
  if (n >= 0) probability = n / (double) size;
  seed_seq seq{DetProb::seed};
  vector<unsigned int> local_seeds (kPCore);
  seq.generate(local_seeds.begin(), local_seeds.end());
  long long global_size = 0;
  vector<string> cols = det_sql.att_cols;
  cols.insert(cols.begin(), det_sql.obj_col);
  int att_count = (int) det_sql.att_cols.size();
  MeanVar mv = MeanVar(att_count);
  int m = det_sql.att_cols.size() + (det_sql.has_count_constraint);
  #pragma omp parallel num_threads(kPCore)
  {
    int seg = omp_get_thread_num();
    default_random_engine gen (local_seeds[seg]);
    uniform_real_distribution dist (0.0, 1.0);
    long long start_id = seg * chunk + 1;
    long long end_id = (seg + 1) * chunk;
    RMatrixXd local_A; vector<long long> local_ids;
    PgManager *_pg = new PgManager();
    if (probability >= 1){
      _pg->getConsecutiveTuples(local_A, local_ids, det_sql.table_name, start_id, end_id, cols, det_sql.filter_cols, det_sql.filter_intervals);
    } else{
      vector<long long> selected_ids; selected_ids.reserve((int)(chunk * probability));
      for (long long id = start_id; id <= end_id; id ++){
        if (dist(gen) <= probability) selected_ids.push_back(id);
      }
      _pg->getSelectedTuples(local_A, local_ids, det_sql.table_name, selected_ids, cols, det_sql.filter_cols, det_sql.filter_intervals);
    }
    delete _pg;
    MeanVar local_mv = MeanVar(att_count);
    int local_size = (int) local_ids.size();
    long long local_start = -1;
    #pragma omp critical (c2)
    {
      local_start = global_size;
      global_size += local_size;
    }
    #pragma omp barrier
    #pragma omp master
    {
      resize(m, global_size);
    }
    #pragma omp barrier
    VectorXd x (att_count);
    for (int i = 0; i < local_size; i ++){
      for (int j = 0; j < att_count; j ++){
        x(j) = local_A(j+1, i);
      }
      local_mv.add(x);
    }
    #pragma omp critical (c3)
    {
      mv.add(local_mv);
    }
    if (local_size > 0){
      memcpy(&(ids[local_start]), &(local_ids[0]), local_size*sizeof(long long));
      memcpy(&c(local_start), &local_A(0, 0), local_size*sizeof(double));
      for (int j = 0; j < att_count; j ++){
        memcpy(&A(j, local_start), &local_A(j+1, 0), local_size*sizeof(double));
      }
    }
  }
  if (!det_sql.is_maximize) c = -c;
  if (det_sql.has_count_constraint) std::fill(&A(att_count, 0), (&A(att_count, 0))+global_size, 1.0);
  det_bound = DetBound(det_sql.att_senses, mv.getMean(), mv.getVar(), DetProb::seed);
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

double DetProb::generateBounds(double E, double alpha, double hardness){
  bl.fill(-DBL_MAX);
  bu.fill(DBL_MAX);
  if (det_sql->has_count_constraint){
    bl(bl.size() - 1) = E*kLowerCountFactor;
    bu(bu.size() - 1) = E*kUpperCountFactor;
  }
  return det_bound.sampleHardness(E, alpha, hardness, bl, bu);
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
  DetProb::seed = seed;
  if (seed < 0) {
    random_device rd;
    DetProb::seed = rd();
  }
  det_bound.setSeed(seed);
}

void DetProb::copyBounds(VectorXd &bl, VectorXd &bu, double cl, double cu){
  assert(DetProb::bl.size() > bl.size());
  assert(DetProb::bu.size() > bu.size());
  int att_count = (int) bl.size();
  if (att_count == 0) return;
  memcpy(&(DetProb::bl(0)), &bl(0), att_count*sizeof(double));
  memcpy(&(DetProb::bu(0)), &bu(0), att_count*sizeof(double));
  DetProb::bl(att_count) = cl;
  DetProb::bu(att_count) = cu;
}

void DetProb::display(){
  cout << "----------------- DetProb Configurations -----------------" << endl;
  cout << "------ u (of each tuple, assume same for all tuples) ------" << endl;
  cout << this->u(0) << endl;
  cout << "------- bl & bu (bound of constraints incl. cl &cu) -------" << endl;
  cout << this->bl << endl;
  cout << this->bu << endl; 
  cout << "-----------------------------------------------------------" << endl;
}