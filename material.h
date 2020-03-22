#ifndef __MATERIAL_H__
#define __MATERIAL_H__

#include "vec3.h"
#include "ray.h"
#include "hitable.h"

// 终于到了材质这一步
class Material {
 public:
  virtual bool Scatter(const Ray& ray_in, const HitRecord& rec,
                       Vec3* attenuation, Ray* scattered) const = 0;
};

#endif  // __MATERIAL_H__
