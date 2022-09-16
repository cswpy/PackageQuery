#pragma once

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

const string kProjectHome = configs["project_home"];

static const int kPCore = stoi(configs["physical_core"]);
static const int kLCore = stoi(configs["logical_core"]);
static const int kMaxMultiColumnIndexes = stoi(configs["max_multicolumn_indexes"]);
static const int kPgPort = stoi(configs["port"]);
static const int kPrecision = stoi(configs["precision"]);
static const long long kMainMemorySize = stoll(configs["main_memory_size"]);
static const long long kInMemorySize = stoll(configs["in_memory_size"]); // In Memory Size for all cores in term of number of tuples

const string kPgDatabase = configs["database"];
const string kPgUser = configs["user"];
const string kPgPassword = configs["password"];
const string kPgHostaddr = configs["hostname"];

const string kId = configs["id_column"]; // Must not equal to "tid"
const string kStatTable = "dlv_stats";
const string kPartitionTable = "dlv_partitions";
const string kNullLiteral = "null";

long long getTupleCount(int tuple_size, int core_count=kPCore);