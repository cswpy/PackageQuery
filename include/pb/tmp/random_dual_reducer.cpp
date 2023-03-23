#include "pb/core/gurobi_solver.h"
#include "pb/core/parallel_pq.h"
#include "pb/util/unumeric.h"
#include "pb/lib/common.h"
#include "random_dual_reducer.h"

#include "pb/util/udebug.h"

static const int kIlpSize = 5000;
static const double kEpsilon = 1e-8;
static const double kSafeMipGap = DBL_MAX;
enum StayStatus {NotStay, LikelyStay, Stay};

RandomDualReducer::~RandomDualReducer(){
}

DetProb* RandomDualReducer::filtering(VectorXi &reduced_index, int core, const DetProb &prob, VectorXd &dual_sol, int stay_count, int stay_mode, vector<int> &stay){
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
      // cout << local_start_index << " " << start_stay_count << endl;
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

RandomDualReducer::RandomDualReducer(int core, const DetProb &prob, bool is_safe, int max_failure, double min_gap, double time_limit){
  failure_count = 0;
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
    // cout << "OK1\n";
    Dual dual = Dual(core, prob);
    if (dual.status == Found){
      memcpy(&(lp_sol(0)), &(dual.sol(0)), n*sizeof(double));
      lp_score = dual.score;
      exe_lp = dual.exe_solve;
      vector<int> stay (n, NotStay);
      int stay_count = 0;
      int stay_mode = LikelyStay;
      for (int i = 0; i < m; i ++){
        if (dual.bhead(i) < n){
          stay[dual.bhead(i)] = Stay;
          stay_count ++;
        }
      }
      for (int i = 0; i < n; i ++){
        if (isGreater(dual.sol(i), prob.l(i))){
          if (stay[i] != Stay){
            stay[i] = Stay;
            stay_count ++;
          }
        }
      }
      vector<int> inds (n);
      iota(inds.begin(), inds.end(), 0);
      std::random_shuffle(inds.begin(), inds.end());
      for (int i = 0; i < n; i ++){
        if (stay[inds[i]] != Stay){
          stay[inds[i]] = Stay;
          stay_count ++;
          if (stay_count == kIlpSize) break;
        }
      }
      // Filtering
      VectorXi reduced_index (stay_count);
      DetProb *reduced_prob = filtering(reduced_index, core, prob, dual.sol, stay_count, stay_mode, stay);
      // cout << "OK3\n";
      GurobiSolver gs = GurobiSolver(*reduced_prob);
      gs.solveIlp(min_gap, 5.0);
      exe_gb += gs.exe_ilp;
      status = gs.ilp_status;
      // cout << "OK4\n";
      if (gs.ilp_status == Found){
        for (int i = 0; i < gs.ilp_sol.size(); i ++){
          // cout << i << " " << gs.ilp_sol.size() << " " << reduced_index.size() << " " << reduced_prob->c.size() << endl;
          // cout << "HERE " << ilp_sol.size() << " " << reduced_index(i) << endl;
          ilp_sol(reduced_index(i)) = gs.ilp_sol(i);
          ilp_score += gs.ilp_sol(i) * reduced_prob->c(i);
        }
        // cout << "OK5\n";
        delete reduced_prob;
      } else if (is_safe && gs.ilp_status != Timeout){
        int current_size = kIlpSize;
        vector<int> inds (n);
        iota(inds.begin(), inds.end(), 0);
        random_shuffle(inds.begin(), inds.end());
        int current_ind = -1;
        while (true){
          failure_count ++;
          if (max_failure >= 0){
            if (failure_count > max_failure){
              delete reduced_prob;
              break;
            }
          }
          delete reduced_prob;
          current_size = min(current_size*2, n);
          while (current_ind < n){
            current_ind ++;
            if (stay[inds[current_ind]] != Stay){
              stay[inds[current_ind]] = Stay;
              stay_count ++;
              if (stay_count == current_size) break;
            }
          }
          // cout << "SAFE " << stay_count << endl;
          reduced_prob = filtering(reduced_index, core, prob, dual.sol, stay_count, stay_mode, stay);
          GurobiSolver _gs = GurobiSolver(*reduced_prob);
          // cout << "HERE " << current_size << " with " << time_limit << endl;
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