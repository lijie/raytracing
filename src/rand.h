#ifndef __RAND_H__
#define __RAND_H__

#include <random>
#include "vec3.h"

// 返回随机数 [0, 1)
inline double Rand() {
  int v = rand();
  return double(v) / double(RAND_MAX);
}

// [start, end)
double RandDouble(double start = 0, double end = 1);
Vec3 RandomCosineDirection();
Vec3 RandomToSphere(double radius, double distance_squared);

#endif  // __RAND_H__