#ifndef __HITABLE_H__
#define __HITABLE_H__

#include <string>
#include "vec3.h"
#include "ray.h"

class Material;
class Hitable;

struct HitRecord {
  double t;         // paramter of ray
  Vec3 p;           // ray at point
  Vec3 normal;      // normal vector of this point
  Hitable* target;  // hitted target, for debug
  Material* mat;
};

class Hitable {
 public:
  virtual bool Hit(const Ray& ray, double t_min, double t_max,
                   HitRecord* rec) const = 0;
  virtual std::string Name() = 0;
};

#endif
