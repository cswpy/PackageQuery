#pragma once

#include <string>

using std::string;

const string kProjectHome = "/home/alm818/PackageQuery";
//const string kProjectHome = "C:/Users/xuana/Desktop/VisualStudioCode/PackageQuery";

constexpr int kPCore = 80;
constexpr int kLCore = 160;
constexpr int kMaxMultiColumnIndexes = 2;
constexpr int kPgPort = 5433;
constexpr int kPrecision = 20;
constexpr long long kInMemorySize = 1000000; // In Memory Size for all cores in term of number of tuples

const string kPgUser = "alm818";
const string kPgPassword = "";
const string kPgHostaddr = "127.0.0.1";

const string kId = "id"; // Must not equal to "tid"
const string kStatTable = "dlv_stats";
const string kPartitionTable = "dlv_partitions";
const string kNullLiteral = "null";
