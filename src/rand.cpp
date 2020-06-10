#include "rand.h"
#include <random>

static std::random_device rd;
static std::mt19937 mt(rd());

void RandInit() {
  // std::random_device rd;
  // std::mt19937 mt(rd());
  // std::uniform_real_distribution<double> dist;
}

double RandDouble(double start, double end) {
    std::uniform_real_distribution<double> dist(start, end);
    return dist(mt);
}

Vec3 RandomCosineDirection() {
  double r1 = RandDouble();
  double r2 = RandDouble();
  double z = sqrt(1 - r2);
  double phi = 2 * M_PI * r1;
  double x = cos(phi) * sqrt(r2);
  double y = sin(phi) * sqrt(r2);
  return Vec3(x, y, z);
}

