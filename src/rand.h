#ifndef __RAND_H__
#define __RAND_H__

#include <random>

// 返回随机数 [0, 1)
inline double Rand() {
  int v = rand();
  return double(v) / double(RAND_MAX);
}

// [start, end)
double RandDouble(double start, double end);
#endif  // __RAND_H__