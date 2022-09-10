#pragma once

#include "pb/util/udebug.h"

class Synthetic{
public:
  static string table_name;
  Profiler pro;
private:
  void init();
  void create(long long N, int ucount, int ncount, vector<double> means, vector<double> vars);
public:
  Synthetic();
  void createUniform(long long N, int count, double mean, double var);
  void createNormal(long long N, int count, double mean, double var);
  void createMixed(long long N, int ucount, int ncount, double mean, double var);
};