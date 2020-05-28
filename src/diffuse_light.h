#ifndef __SRC_DIFFUSE_LIGHT_H__
#define __SRC_DIFFUSE_LIGHT_H__

#include "material.h"
#include "texture.h"

class DiffuseLight : public Material {
 public:
  DiffuseLight(Texture* tex) : tex_(tex) {}
  bool Scatter(const Ray& ray_in, const HitRecord& rec, Vec3* attenuation,
               Ray* scattered) const override {
    return false;
  }

  Vec3 Emmited(double u, double v, const Vec3& point) const {
    return tex_->Value(u, v, point);
  }

 private:
  Texture* tex_;
};

#endif  // __SRC_DIFFUSE_LIGHT_H__