#ifndef __ORTHONORMAL_BASIS__
#define __ORTHONORMAL_BASIS__

#include "vec3.h"

class OrthonormalBasis {
 public:
  OrthonormalBasis() {}
  Vec3 u() const { return axis_[0]; }
  Vec3 v() const { return axis_[1]; }
  Vec3 w() const { return axis_[2]; }

  Vec3 Local(const Vec3& a) const { return a.x() * u() + a.y() * v() + a.z() * w(); }
  Vec3 Local(double x, double y, double z) const {
    return x * u() + y * v() + z * w();
  }

  void BuildFromW(const Vec3& n) {
    axis_[2] = unit_vector(n);
    Vec3 a = fabs(w().x()) > 0.9 ? Vec3(0, 1, 0) : Vec3(1, 0, 0);
    axis_[1] = unit_vector(cross(w(), a));
    axis_[0] = unit_vector(cross(v(), w()));
  }

  // see pbrt-book ch2.2.4  coordinate system from a vector
  void BuildFromW2(const Vec3& n) {}

 private:
  Vec3 axis_[3];
};

#endif  // __ORTHONORMAL_BASIS__
