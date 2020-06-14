#ifndef __MATERIAL_H__
#define __MATERIAL_H__

#include <memory>

#include "hitable.h"
#include "ray.h"
#include "vec3.h"

class Pdf;

struct ScatterRecord {
  Ray specular_ray;
  bool is_specular;
  Vec3 attenuation;
  std::shared_ptr<Pdf> pdf;
};

// 终于到了材质这一步
class Material {
 public:
  virtual bool Scatter(const Ray& ray_in, const HitRecord& rec,
                       Vec3* attenuation, Ray* scattered,
                       double* pdf) const = 0;

  virtual bool Scatter(const Ray& ray_in, const HitRecord& rec,
                       ScatterRecord* srec) const {
    return false;
  }

  virtual double ScatteringPdf(const Ray& ray, const HitRecord& rec,
                               const Ray& scattered) const {
    return 1.0;
  }
  // 自发光
  virtual Vec3 Emmited(double u, double v, const Vec3& point) const {
    return Vec3(0, 0, 0);
  }
};

#endif  // __MATERIAL_H__
