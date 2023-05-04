#include "pb/core/gurobi_solver.h"
#include "pb/core/parallel_pq.h"
#include "pb/util/unumeric.h"
#include "test_dual_reducer.h"

#include "pb/util/udebug.h"

static const int kIlpSize = 5000;
static const double kEpsilon = 1e-8;
static const double kSafeMipGap = DBL_MAX;
enum StayStatus {NotStay, LikelyStay, Stay};

TestDualReducer::~TestDualReducer(){
}

DetProb* TestDualReducer::filtering(VectorXi &reduced_index, int core, const DetProb &prob, VectorXd &dual_sol, int stay_count, int stay_mode, vector<int> &stay){
  int n = (int) prob.c.size();
  int m = (int) prob.bl.size();  
  reduced_index.resize(stay_count);
  DetProb *reduced_prob = new DetProb(m, stay_count);
  reduced_prob->bl = prob.bl;
  reduced_prob->bu = prob.bu;
  int start_stay_count = 0;
  #pragma omp parallel num_threads(core)
  {
    int local_start_index = -1;
    int local_stay_count = 0;
    VectorXi local_stay_index (stay_count);
    VectorXd local_change_b (m); local_change_b.fill(0);
    #pragma omp for nowait
    for (int i = 0; i < n; i ++){
      if (stay[i] == Stay || stay[i] == stay_mode){
        local_stay_index(local_stay_count) = i;
        local_stay_count ++;
      } else{
        ilp_sol(i) = dual_sol(i);
        #pragma omp atomic
        ilp_score += dual_sol(i) * prob.c(i);
        for (int j = 0; j < m; j ++){
          local_change_b(j) += prob.A(j, i) * dual_sol(i);
        }
      }
    }
    #pragma omp critical (c3)
    {
      local_start_index = start_stay_count;
      start_stay_count += local_stay_count;
    }
    if (local_stay_count) memcpy(&(reduced_index(local_start_index)), &(local_stay_index(0)), local_stay_count*sizeof(int));
    // for (int i = local_start_index; i < local_start_index + local_stay_count; i ++){
    //   reduced_index(i) = local_stay_index(i-local_start_index);
    //   assert(reduced_index(i) == local_stay_index(i-local_start_index));
    // }
    #pragma omp barrier
    for (int i = 0; i < m; i ++){
      if (prob.bl(i) != -DBL_MAX){
        #pragma omp atomic
        reduced_prob->bl(i) -= local_change_b(i);
      }
      if (prob.bu(i) != DBL_MAX){
        #pragma omp atomic
        reduced_prob->bu(i) -= local_change_b(i);
      }
    }
  }
  for (int i = 0; i < start_stay_count; i ++){
    for (int j = 0; j < m; j ++){
      reduced_prob->A(j, i) = prob.A(j, reduced_index(i));
    }
    reduced_prob->c(i) = prob.c(reduced_index(i));
    reduced_prob->l(i) = prob.l(reduced_index(i));
    reduced_prob->u(i) = prob.u(reduced_index(i));
  }
  return reduced_prob;
}

TestDualReducer::TestDualReducer(int core, const DetProb &prob, bool is_safe, double min_gap, double time_limit){
  std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
  int n = (int) prob.c.size();
  int m = (int) prob.bl.size();
  ilp_sol.resize(n); lp_sol.resize(n);
  ilp_score = lp_score = exe_ilp = exe_lp = exe_gb = 0;
  status = NotFound;
  if (n < kIlpSize){
    // Gurobi
    GurobiSolver gs = GurobiSolver(prob);
    gs.solveLp();
    exe_lp = gs.exe_lp;
    lp_sol = gs.lp_sol;
    lp_score = gs.lp_score;
    gs.solveIlp(min_gap, time_limit);
    exe_gb += gs.exe_ilp;
    status = gs.ilp_status;
    if (gs.ilp_status == Found){
      ilp_sol = gs.ilp_sol;
      ilp_score = gs.ilp_score;
    }
  } else{
    Dual dual = Dual(core, prob);
    if (dual.status == Found){
      memcpy(&(lp_sol(0)), &(dual.sol(0)), n*sizeof(double));
      lp_score = dual.score;
      exe_lp = dual.exe_solve;
      vector<int> stay (n, NotStay);
      int stay_count = 0;
      int stay_mode = LikelyStay;
      VectorXd left_space (m); left_space.fill(-1);
      VectorXd right_space (m); right_space.fill(-1);
      VectorXd score (n); score.fill(DBL_MAX);
      VectorXd global_a (m); global_a.fill(0);
      int excess_count = 0;
      double max_abs_c = 0;
      vector<pair<double, int>>* excess;
      vector<pair<double, int>>* deficit = new vector<pair<double, int>>(n); // Destruction is handled by ParallelPQ
      for (int i = 0; i < m; i ++){
        if (dual.bhead(i) < n){
          stay[dual.bhead(i)] = Stay;
          stay_count ++;
        }
      }
      int bhead_count = stay_count;
      #pragma omp parallel num_threads(core)
      {
        VectorXd local_a (m); local_a.fill(0);
        double local_max_abs_c = 0;
        #pragma omp for nowait
        for (int j = 0; j < n; j ++){
          if (isEqual(dual.sol(j), prob.u(j), kEpsilon)){
            if (stay[j] == NotStay){
              stay[j] = LikelyStay;
              #pragma omp atomic
              stay_count ++;
            }
            for (int i = 0; i < m; i ++){
              local_a(i) += prob.A(i, j) * prob.u(j);
            }
          }
          local_max_abs_c = max(local_max_abs_c, fabs(prob.c(j)));
        }
        #pragma omp critical (c1)
        {
          max_abs_c = max(max_abs_c, local_max_abs_c);
        }
        for (int i = 0; i < m; i ++){
          #pragma omp atomic
          global_a(i) += local_a(i);
        }
        #pragma omp barrier
        #pragma omp master
        {
          if (stay_count > kIlpSize / 2){
            excess = new vector<pair<double, int>>(stay_count - bhead_count);
            stay_mode = Stay;
          }
          for (int i = 0; i < m; i ++){
            if (prob.bl(i) != -DBL_MAX) left_space(i) = global_a(i) - prob.bl(i);
            if (prob.bu(i) != DBL_MAX) right_space(i) = prob.bu(i) - global_a(i);
          }
        }
        #pragma omp barrier
        #pragma omp for nowait
        for (int i = 0; i < n; i ++){
          if (isLess(dual.sol(i), prob.u(i), kEpsilon)){
            score(i) = 0;
            for (int j = 0; j < m; j ++){
              if (prob.A(j, i) > 0){
                if (isGreater(right_space(j), 0, kEpsilon)) score(i) += prob.A(j, i) / right_space(j);
                else if (isLess(right_space(j), 0, kEpsilon)) score(i) += (right_space(j) - prob.A(j, i)) / right_space(j);
              } else if (prob.A(j, i) < 0){
                if (isGreater(left_space(j), 0, kEpsilon)) score(i) -= prob.A(j, i) / left_space(j);
                else if (isLess(left_space(j), 0, kEpsilon)) score(i) += (left_space(j) + prob.A(j, i)) / left_space(j);
              }
            }
            // if (m > 0) score(i) /= m;
            if (max_abs_c > 0) score(i) -= prob.c(i) / max_abs_c;
          }
          (*deficit)[i] = {-score(i), i};
        }

        if (stay_count > kIlpSize / 2){
          vector<int> local (stay_count);
          int local_count = 0;
          #pragma omp for nowait
          for (int i = 0; i < n; i ++){
            if (stay[i] == LikelyStay){
              local[local_count] = i;
              local_count ++;
            }
          }
          int local_start = -1;
          #pragma omp critical (c2)
          {
            local_start = excess_count;
            excess_count += local_count;
          }
          for (int i = 0; i < local_count; i ++){
            (*excess)[local_start + i] = {-prob.c(local[i]), local[i]};
          }
        }
      }
      if (stay_count > kIlpSize / 2){
        int excess_pop = kIlpSize / 2 - bhead_count;
        int zero_count = n - stay_count;
        if (kIlpSize / 2 + zero_count < kIlpSize){
          excess_pop = kIlpSize - zero_count - bhead_count;
        }
        ParallelPQ excess_pq = ParallelPQ(core, excess);
        for (int i = 0; i < excess_pop; i ++){
          stay[excess_pq.peak().second] = Stay;
          excess_pq.pop();
        }
        stay_count = excess_pop + bhead_count;
      }
      ParallelPQ deficit_pq = ParallelPQ(core, deficit);
      while (stay_count < kIlpSize){
        int index = deficit_pq.peak().second;
        if (stay[index] != Stay){
          stay[index] = Stay;
          stay_count ++;
        }
        deficit_pq.pop();
      }
      // Filtering
      VectorXi reduced_index (stay_count);
      DetProb *reduced_prob = filtering(reduced_index, core, prob, dual.sol, stay_count, stay_mode, stay);
      GurobiSolver gs = GurobiSolver(*reduced_prob);
      gs.solveIlp(min_gap, 5.0);
      exe_gb += gs.exe_ilp;
      status = gs.ilp_status;
      if (gs.ilp_status == Found){
        for (int i = 0; i < gs.ilp_sol.size(); i ++){
          ilp_sol(reduced_index(i)) = gs.ilp_sol(i);
          ilp_score += gs.ilp_sol(i) * reduced_prob->c(i);
        }
        delete reduced_prob;
      } else if (is_safe && gs.ilp_status != Timeout){
        int current_size = kIlpSize;
        while (true){
          delete reduced_prob;
          current_size = min(current_size*2, n);
          // cout << "AT " << current_size << endl;
          while (stay_count < current_size){
            int index = deficit_pq.peak().second;
            if (stay[index] != Stay){
              stay[index] = Stay;
              stay_count ++;
            }
            deficit_pq.pop();
          }
          reduced_prob = filtering(reduced_index, core, prob, dual.sol, stay_count, stay_mode, stay);
          GurobiSolver _gs = GurobiSolver(*reduced_prob);
          _gs.solveIlp(kSafeMipGap, time_limit);
          exe_gb += _gs.exe_ilp;
          status = _gs.ilp_status;
          if (_gs.ilp_status == Found){
            for (int i = 0; i < _gs.ilp_sol.size(); i ++){
              ilp_sol(reduced_index(i)) = _gs.ilp_sol(i);
              ilp_score += _gs.ilp_sol(i) * reduced_prob->c(i);
            }
            delete reduced_prob;
            break;
          } else if (current_size == n || _gs.ilp_status == Timeout){
            delete reduced_prob;
            break;
          }
        }
      } else{
        // Not safe, LP found but no ILP found 
        status = HalfFeasible;
        ilp_sol = lp_sol;
        ilp_score = lp_score;
      }
    } else{
      status = dual.status;
    }
  }
  exe_ilp = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start).count() / 1000000.0;
}