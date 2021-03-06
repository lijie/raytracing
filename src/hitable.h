#ifndef __HITABLE_H__
#define __HITABLE_H__

#include <string>

#include "ray.h"
#include "vec3.h"

class Material;
class Hitable;
class AABB;
class Pdf;

struct HitRecord {
  double t;         // paramter of ray
  Vec3 p;           // ray at point
  Vec3 normal;      // normal vector of this point
  double u;
  double v;
  Hitable* target;  // hitted target, for debug
  Material* mat;
};

class Hitable {
 public:
  virtual bool Hit(const Ray& ray, double t_min, double t_max,
                   HitRecord* rec) const = 0;
  virtual bool BoundingBox(double t0, double t1, AABB* box) const = 0;
  virtual Vec3 SurfaceRandomPoint(const Vec3& origin) {
    return Vec3(0, 0, 0);
  }
  virtual double SurfaceRandomPdf(const Vec3& origin, const Vec3& direction) {
    return 0;
  }
  virtual std::string Name() = 0;
};

#endif
