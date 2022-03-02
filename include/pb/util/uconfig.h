#pragma once

#include <string>

using std::string;

const string kProjectHome = "/home/alm818/PackageQuery";
//const string kProjectHome = "C:/Users/xuana/Desktop/VisualStudioCode/PackageQuery";

constexpr int kPCore = 80;
constexpr int kLCore = 160;
constexpr int kMaxMultiColumnIndexes = 32;
constexpr int kPgPort = 5433;
constexpr int kPrecision = 20;
constexpr int kInMemorySize = 1000000;

const string kPgUser = "alm818";
const string kPgPassword = "";
const string kPgHostaddr = "127.0.0.1";

const string kId = "id";
const string kStatTable = "dlv_stats";
const string kPartitionTable = "dlv_partitions";
const string kNullLiteral = "null";
