#pragma once

constexpr double kNumericEps = 1e-12;

bool isEqual(double x, double y, double eps=kNumericEps);
bool isLess(double x, double y, double eps=kNumericEps);
bool isGreater(double x, double y, double eps=kNumericEps);
bool isLessEqual(double x, double y, double eps=kNumericEps);
bool isGreaterEqual(double x, double y, double eps=kNumericEps);
bool isInteger(double x, double eps=kNumericEps);

long long ceilDiv(long long x, long long q);

template <typename T>
double sign(T value) {
  if (value < 0) return -1.0;
  if (value > 0) return 1.0;
  return 0.0;
}

template <typename T>
double nonNegativeSign(T value) {
  if (value < 0) return -1.0;
  return 1.0;
}