// axis-aligned bouding box
#ifndef __AABB_H__
#define __AABB_H__

#include "ray.h"
#include "vec3.h"

#define fmin(a, b) (a) <= (b) ? (a) : (b)
#define fmax(a, b) (a) >= (b) ? (a) : (b)

class AABB {
 public:
  AABB() {}
  AABB(const Vec3& min, const Vec3& max) : min_(min), max_(max) {}

  Vec3 min() const { return min_; }
  Vec3 max() const { return max_; }

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

  static AABB SurroudingBox(const AABB& box1, const AABB& box2) {
    Vec3 min(fmin(box1.min().x(), box2.min().x()),
             fmin(box1.min().y(), box2.min().y()),
             fmin(box1.min().z(), box2.min().z()));
    Vec3 max(fmax(box1.max().x(), box2.max().x()),
             fmax(box1.max().y(), box2.max().y()),
             fmax(box1.max().z(), box2.max().z()));
    return AABB(min, max);
  }

 private:
  Vec3 min_;
  Vec3 max_;
};
#endif
