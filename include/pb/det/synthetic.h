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
  void create(string method, long long N, vector<double> args);
public:
  Synthetic(string dbname);
  void createUniform(long long N, double mean, double var);
  void createNormal(long long N, double mean, double var);
};