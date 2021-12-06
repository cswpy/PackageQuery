#pragma once

#include "pb/util/udebug.h"

class DynamicLowVariance{
private:
  string dbname;
public:
  static double kGroupRatio;
  static int kMinPartitionSize;
  static int kInitialSize;

  Profiler pro;
private:
  void init();
public:
  DynamicLowVariance(string dbname);
  void partition();
  void computeStats(string table_name);
};