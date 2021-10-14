#include <iostream>
#include "coo_sparse.h"

CooSparse::CooSparse(){
}

CooSparse::CooSparse(int n, int m){
  this->n = n;
  rows.resize(n);
  for (int i = 0; i < n; i ++) rows[i].reserve(m);
}

void CooSparse::addEntry(int r, int c, double v){
  rows[r].emplace_back(c, v);
}

int CooSparse::getSize(){
  int size = 0;
  for (int i = 0; i < n; i ++) size += rows[i].size();
  return size;
}

VectorXd CooSparse::vectorProduct(VectorXd x){
  VectorXd res (n); res.fill(0);
  for (int i = 0; i < n; i ++){
    for (const auto& pi : rows[i]){
      res[i] += x(pi.first) * pi.second;
    }
  }
  return res;
}