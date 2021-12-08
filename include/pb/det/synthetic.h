#pragma once

#include "pb/util/udebug.h"

class Synthetic{
private:
  string dbname;
public:
  static string table_name;
  static int in_memory_chunk;
  Profiler pro;
private:
  void init();
  void create(long long N, int ucount, int ncount, vector<double> means, vector<double> vars);
public:
  Synthetic(string dbname);
  void createUniform(long long N, int count, double mean, double var);
  void createNormal(long long N, int count, double mean, double var);
  void createMixed(long long N, int ucount, int ncount, double mean, double var);
};