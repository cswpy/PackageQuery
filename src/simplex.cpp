#define FMT_HEADER_ONLY

#include <float.h>
#include <iostream>
#include "utility.h"
#include "simplex.h"
#include "fmt/core.h"

using namespace std;

#define at(row, col) (numcols*(row)+(col))
#define slack(index) (n+(index))
#define artif(index) (npm+1+(index))

void draw(double* tableau, int m, int n){
  cout << "--------------------------------" << endl;
  for (int i = 0; i < m; i ++){
    for (int j = 0; j < n; j ++){
      fmt::print("{: 8.3Lf} ", tableau[i*n+j]);
    }
    cout << endl;
  }
  cout << "--------------------------------" << endl;
}

Simplex::~Simplex(){
  delete[] bhead;
  delete[] tableau;
}

void Simplex::pivot(int horizon, const VectorXd& u, int ent_col){
  #pragma omp master
  {
    lea_row = -1;
    min_ratio = DBL_MAX;
    double col_obj = tableau[at(0, ent_col)];
    if (col_obj < 0){
      // Increase entering variable
      // Check upper bound on entering variable
      if (ent_col < n && u(ent_col) >= 0){
        double ratio = u(ent_col) - tableau[at(1, ent_col)];
        if (min_ratio > ratio){
          min_ratio = ratio;
          lea_row = 1;
        }
      }
    } else{
      // Decrease entering variable
      // Check lower bound on entering variable
      double ratio = tableau[at(1, ent_col)];
      if (min_ratio > ratio){
        min_ratio = ratio;
        lea_row = 1;
      }
    }
    // cout << min_ratio << " HERE" << endl;
    for (int j = 2; j <= m+1; j ++){
      int basic_col = bhead[j-2];
      //cout << "HERE " << j << " " << basic_col << endl;
      double step = sign(col_obj)*tableau[at(j, ent_col)];
      double cur = tableau[at(1, basic_col)];
      // cout << "STEP " << step << " "  << cur << " JENT " << tableau[at(j, ent_col)] << " " << sign(col_obj) << endl;
      double ratio = -1;
      if (basic_col < n){
        if (isGreater(step, 0) && u(basic_col) >= 0) ratio = (u(basic_col) - cur) / step;
        else if (isLess(step,0)) ratio = -cur / step;
      } else{
        if (isLess(step,0)) ratio = -cur / step;
      }
      // if (ratio > 0){
      //   cout << j << " " << ratio << endl;
      // } else if (ratio == 0){
      //   cout << "CAUTION " << j << " " << ratio << endl;
      // }
      if (isGreaterEqual(ratio,0) && min_ratio > ratio){
        min_ratio = ratio;
        lea_row = j;
      }
    }
  }
  #pragma omp barrier
  // #pragma omp single
  // cout << "LEA " << lea_row << " " << min_ratio << " " << endl;
  if (lea_row == 1){
    #pragma omp master
    {
      // cout << "HERE " << ent_col << endl;
      double col_obj = tableau[at(0, ent_col)];
      double step = min_ratio * sign(col_obj);
      tableau[at(1, ent_col)] -= step;
      // No pivot, just change non-basic from lb to ub or ub to lb
      for (int j = 2; j <= m+1; j ++){
        // Update cur
        int basic_col = bhead[j-2];
        double tmp = step*tableau[at(j, ent_col)];
        tableau[at(1, basic_col)] += tmp;
        // if (abs(tableau[at(1, basic_col)]) < kFloatEps) tableau[at(1, basic_col)] = 0;
        // if (tableau[at(1, basic_col)] < 0){
        //   cout << "REALLY BAD "<< tableau[at(1, basic_col)] << " " << j << " " << basic_col << " " << lea_row << " " << min_ratio << " " << endl;
        // }
        // Update rhs
        // tableau[at(j, npm)] += tmp;
      }
    }
  } else if (lea_row == -1){
    #pragma omp master
    {
      status = LS_UNBOUNDED;
    }
  } else{
    // Pivot happening, entering variable ent_col, leaving variable lea_row
    #pragma omp master
    {
      double col_obj = tableau[at(0, ent_col)];
      double step = min_ratio * sign(col_obj);
      tableau[at(1, ent_col)] -= step;
      for (int j = 2; j <= m+1; j ++){
        // Update cur
        int basic_col = bhead[j-2];
        tableau[at(1, basic_col)] += step*tableau[at(j, ent_col)];
        // if (abs(tableau[at(1, basic_col)]) < kFloatEps) tableau[at(1, basic_col)] = 0;
        // if (tableau[at(1, basic_col)] < 0){
        //   cout << "REALLY BAD 2| " << tableau[at(j, ent_col)] << " " << step << " HOLD " << (tableau[at(1, basic_col)]-step*tableau[at(j, ent_col)]) << "+" << step*tableau[at(j, ent_col)] << "=" << tableau[at(1, basic_col)] << " " << j << " " << basic_col << " " << lea_row << " " << min_ratio << " " << endl;
        // }
      }
      bhead[lea_row-2] = ent_col;
      d_cache = tableau[at(lea_row, ent_col)];
      v_cache.resize(m+2);
      for (int j = 0; j <= m+1; j ++) v_cache[j] = tableau[at(j, ent_col)];
    }
    #pragma omp barrier
    #pragma omp for
    for (int i = 0; i < horizon; i ++){
      tableau[at(lea_row, i)] /= d_cache;
    }
    #pragma omp for nowait
    for (int i = 0; i < horizon; i ++){
      for (int j = m+1; j >= 0; j --){
        if (j != 1 && j != lea_row) 
          tableau[at(j, i)] -= tableau[at(lea_row, i)] * v_cache[j];
      }
    }
    // #pragma omp for
    // for (int i = 0; i < horizon; i ++){
    //   if (i != ent_col){
    //     for (int j = m+1; j >= 0; j --){
    //       if (j != 1 && j != lea_row) 
    //         tableau[at(j, i)] -= tableau[at(lea_row, i)] * tableau[at(j, ent_col)];
    //     }
    //   }
    // }
    // #pragma omp master
    // {
    //   for (int j = m+1; j >= 0; j --){
    //     if (j != 1 && j != lea_row) tableau[at(j, ent_col)] = 0;
    //   }
    // }
  }
  #pragma omp barrier
  // #pragma omp single
  // draw(tableau, m+2, numcols);
  //cout << ent_col << endl;
  //cout << min_ratio << " " << lea_row << endl;
}

void Simplex::selectEnteringColumn(int horizon, const VectorXd& u, double& ent_value, int& ent_col){
  // Choosing entering variable
  #pragma omp master
  {
    ent_value = 0;
    ent_col = -1;
  }
  #pragma omp barrier
  double local_ent_value = 0;
  int local_ent_col = -1;
  #pragma omp for nowait
  for (int i = 0; i < horizon; i ++){
    if (i == npm) continue;
    double cur = tableau[at(1, i)];
    double c = tableau[at(0, i)];
    // If the entering variable already lower bounded and obj is non negative
    if (isEqual(cur, 0) && isGreaterEqual(c, 0)) continue;
    // If the entering variable already upper bounded and obj is non positive
    if (i < n && isEqual(cur, u(i)) && isLessEqual(c, 0)) continue;
    if (local_ent_value < fabs(c)){
      local_ent_value = fabs(c);
      local_ent_col = i;
    }
  }
  #pragma omp critical
  {
    if (ent_value < local_ent_value){
      ent_value = local_ent_value;
      ent_col = local_ent_col;
    }
  }
  #pragma omp barrier
}

Simplex::Simplex(int core, const MatrixXd& A, const VectorXd& b, const VectorXd& c, const VectorXd& u){
  status = LS_NOT_FOUND;
  n = c.size();
  m = b.size();
  npm = n+m;
  bhead = new double[m];
  vector<int> ma_map;
  ma = 0;
  for (int i = 0; i < m; i ++){
    if (b(i) >= 0){
      bhead[i] = slack(i);
    } else{
      bhead[i] = artif(ma);
      ma_map.push_back(i);
      ma ++;
    }
  }
  numcols = npm+1+ma;
  tableau = new double[(m+2)*numcols];
  double ent_value = 0;
  int ent_col = -1;

  #pragma omp parallel num_threads(core)
  {
    // Phase-I: Initialization
    #pragma omp for
    for (int i = 0; i < numcols; i ++){
      for (int j = m+1; j >= 0; j --){
        if (j > 1){
          int index = j - 2;
          double s = nonNegativeSign(b(index));
          if (i < n) tableau[at(j, i)] = A(index, i) * s;
          else if (i < npm){
            if (i-n == index) tableau[at(j, i)] = s;
            else tableau[at(j, i)] = 0;
          } else if (i == npm) tableau[at(j, i)] = s * b(index);
          else{
            int ma_index = ma_map[i-npm-1];
            if (ma_index == index) tableau[at(j, i)] = 1;
            else tableau[at(j, i)] = 0;
          }
        } else{
          if (j == 1){
            if (i >= n && i < npm && b(i-n) >= 0) tableau[at(j, i)] = b(i-n);
            else if (i > npm) tableau[at(j, i)] = -b(ma_map[i-npm-1]);
            else tableau[at(j, i)] = 0;
          } else{
            if (i <= npm){
              double val = 0;
              for (auto ma_index : ma_map) val -= tableau[at(ma_index+2, i)];
              tableau[at(j, i)] = val;
            } else tableau[at(j, i)] = 0;
          }
        }
        //#pragma omp critical
        //cout << j << " " << i << "=" << tableau[at(j, i)] << " " << omp_get_thread_num() << endl;
      }
    }
    // #pragma omp single
    // draw(tableau, m+2, numcols);
    // #pragma omp single
    // draw(tableau, m+2, numcols);
    // for (int p = 0; p < 25; p ++){
    //   selectEnteringColumn(numcols, u, ent_value, ent_col);
    //   #pragma omp single
    //   cout << p << " ENT " << ent_col << " " << ent_value << endl;
    //   pivot(numcols, u, ent_col);
    // }
    // Phase-I
    while (true){
      selectEnteringColumn(numcols, u, ent_value, ent_col);
      #pragma omp master
      {
        //cout << "ENT " << ent_value << " " << ent_col << endl;
        if (ent_value == 0){
          // cout << "WTF " << tableau[at(0, npm)] << " " << endl;
          if (isGreaterEqual(tableau[at(0, npm)], 0)) status = LS_FEASIBLE;
          else status = LS_INFEASIBLE;
        }
      }
      #pragma omp barrier
      if (status != LS_NOT_FOUND) break;
      pivot(numcols, u, ent_col);
      if (status == LS_UNBOUNDED) break;
    }

    if (status == LS_FEASIBLE){
      // #pragma omp single
      // draw(tableau, m+2, numcols);
      // Phase-II: Initialization
      #pragma omp master
      {
        sum = 0;
      }
      #pragma omp barrier
      double local_sum = 0;
      #pragma omp for nowait
      for (int i = 0; i <= npm; i ++){
        if (i < n) tableau[at(0, i)] = -c(i);
        else tableau[at(0, i)] = 0;
        local_sum += tableau[at(0, i)] * tableau[at(1, i)];
      }
      #pragma omp atomic
      sum += local_sum;
      #pragma omp barrier
      #pragma omp master
      {
        tableau[at(0, npm)] = sum;
        v_cache.resize(m);
        for (int j = 2; j <= m+1; j ++){
          int basic_col = bhead[j-2];
          v_cache[j-2] = tableau[at(0, basic_col)];
        }
      }
      #pragma omp barrier
      // #pragma omp single
      // draw(tableau, m+2, numcols);
      #pragma omp for
      for (int i = 0; i <= npm; i ++){
        double sub_obj = 0;
        for (int j = 2; j <= m+1; j ++){
          int basic_col = bhead[j-2];
          if (basic_col < n) sub_obj += tableau[at(j, i)] * v_cache[j-2];
        }
        tableau[at(0, i)] -= sub_obj;
      }

      // #pragma omp single
      // draw(tableau, m+2, numcols);
      // #pragma omp single
      // {
      //   for (int i = 0; i < m; i ++){
      //     cout << bhead[i] << " ";
      //   }
      //   cout << endl;
      // }
      // #pragma omp single
      // draw(tableau, m+2, numcols);

      // Phase-II
      while (true){
        selectEnteringColumn(npm+1, u, ent_value, ent_col);
        #pragma omp master
        {
          //cout << "ENT " << ent_value << " " << ent_col << endl;
          if (ent_value == 0){
            status = LS_FOUND;
            // cout << "WTF " << tableau[at(0, npm)] << " " << endl;
          }
        }
        #pragma omp barrier
        if (status == LS_FOUND) break;
        pivot(npm+1, u, ent_col);
        if (status == LS_UNBOUNDED) break;

        // #pragma omp single
        // draw(tableau, m+2, numcols);
        // #pragma omp single
        // {
        //   for (int i = 0; i < m; i ++){
        //     cout << bhead[i] << " ";
        //   }
        //   cout << endl;
        // }
      }

      // #pragma omp single
      // draw(tableau, m+2, numcols);
      // #pragma omp single
      // {
      //   for (int i = 0; i < m; i ++){
      //     cout << bhead[i] << " ";
      //   }
      //   cout << endl;
      // }
    }
  }
}