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
