#pragma once

#define DEBUG 1

#include <string>
#include <unordered_map>
#include <unistd.h>
#include <fstream>
#include <iostream>

#include "umisc.h"

using std::string;
using std::unordered_map;

static int readConfigs(unordered_map<string, string>& configs){
  char cwd[PATH_MAX];
  if (getcwd(cwd, sizeof(cwd)) == NULL) {
    perror("getcwd() error");
    return 1;
  }
  string project_home (cwd);
  project_home = project_home.substr(0, project_home.find_last_of(separator()));
  configs["project_home"] = project_home;

  std::fstream config_file; config_file.open(project_home + separator() + "config.txt", std::ios::in);
  if (config_file.is_open()){
    string line;
    while (getline(config_file, line)){
      line = trim(line);
      if (line.length() > 0 && line[0] != '#'){
        auto pos = line.find('=');
        if (pos != string::npos){
          string key = line.substr(0, pos);
          string value = line.substr(pos+1);
          configs[trim(key)] = trim(value);
        }
      }
    }
    config_file.close();
  }
  // for (auto p : configs){
  //   std::cout << p.first << " " << p.second << "\n";
  // }
  return 0;
}

static unordered_map<string, string> configs;
static int _ = readConfigs(configs);

#if DEBUG
  static string kProjectHome = configs["project_home"];

  static int kPCore = stoi(configs["physical_core"]);
  static int kLCore = stoi(configs["logical_core"]);
  static int kMaxMultiColumnIndexes = stoi(configs["max_multicolumn_indexes"]);
  static int kPgPort = stoi(configs["port"]);
  static int kPrecision = stoi(configs["precision"]);
  static int kGlobalSeed = stoi(configs["global_seed"]);
  static long long kLpSize = stoll(configs["lp_size"]);
  static double kMainMemorySize = stod(configs["main_memory_size"]);

  static string kPgDatabase = configs["database"];
  static string kPgUser = configs["user"];
  static string kPgPassword = configs["password"];
  static string kPgHostaddr = configs["hostname"];
  static string kPgSchema = configs["schema"];

  static string kId = configs["id_column"]; // Must not equal to "tid"
  static string kStatTable = "stats";
  static string kPartitionTable = "dlv_partitions";
  static string kNullLiteral = "null";
  static string kIntervalType = "floatrange";
  static string kTempPrefix = "tmp";
#else
  static const string kProjectHome = configs["project_home"];

  static const int kPCore = stoi(configs["physical_core"]);
  static const int kLCore = stoi(configs["logical_core"]);
  static const int kMaxMultiColumnIndexes = stoi(configs["max_multicolumn_indexes"]);
  static const int kPgPort = stoi(configs["port"]);
  static const int kPrecision = stoi(configs["precision"]);
  static const int kGlobalSeed = stoi(configs["global_seed"]);
  static const long long kLpSize = stoll(configs["lp_size"]);
  static const double kMainMemorySize = stod(configs["main_memory_size"]);

  static const string kPgDatabase = configs["database"];
  static const string kPgUser = configs["user"];
  static const string kPgPassword = configs["password"];
  static const string kPgHostaddr = configs["hostname"];
  static const string kPgSchema = configs["schema"];

  static const string kId = configs["id_column"]; // Must not equal to "tid"
  static const string kStatTable = "stats";
  static const string kPartitionTable = "dlv_partitions";
  static const string kNullLiteral = "null";
  static const string kIntervalType = "floatrange";
  static const string kTempPrefix = "tmp";
#endif

long long getTupleCount(int tuple_size, int core_count=kPCore, double main_memory=kMainMemorySize);