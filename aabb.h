// axis-aligned bouding box

#include "ray.h"
#include "vec3.h"

class AABB {
 public:
  AABB(const Vec3& min, const Vec3& max) : min_(min), max_(max) {}

  Vec3 min() { return min_; }
  Vec3 max() { return max_; }

  bool Hit(const Ray& ray, double tmin, double tmax) const {
    for (int i = 0; i < 3; i++) {
      double t0 = (min_[i] - ray.origin()[i]) / ray.direction()[i];
      double t1 = (max_[i] - ray.origin()[i]) / ray.direction()[i];
      if (t0 > t1) std::swap(t0, t1);
      tmin = t0 < tmin ? tmin : t0;
      tmax = t1 > tmax ? tmax : t1;
      if (tmax < tmin) return false;
    }
    return true;
  }

 private:
  Vec3 min_;
  Vec3 max_;
};