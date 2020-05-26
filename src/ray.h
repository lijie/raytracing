#ifndef __RAY_H__
#define __RAY_H__

#include <string>
#include <sstream>
#include "vec3.h"

class Ray {
 public:
  Ray() {}
  Ray(const Vec3& a, const Vec3& b, double time = 0) {
    a_ = a;
    b_ = b;
    time_ = time;
  }

  std::string ToString() const {
    std::stringstream ss;
    ss << a_.ToString();
    ss << " -> ";
    ss << b_.ToString();
    return ss.str();
  }

  Vec3 origin() const { return a_; }
  Vec3 direction() const { return b_; }
  Vec3 point_at_parameter(double t) const { return a_ + t * b_; }
  double time() const { return time_; }
  Vec3 a_;  // 原点
  Vec3 b_;  // 方向
  double time_;
};

#endif
