#ifndef __LAMBERTIAN_H__
#define __LAMBERTIAN_H__

#include "texture.h"
#include "material.h"
#include "rand.h"

class Lambertian : public Material {
 public:
  Lambertian(const Vec3& a) {
    albedo_ = new ConstantTexture(a);
  }
  Lambertian(Texture *tex) : albedo_(tex) {}
  bool Scatter(const Ray& ray_in, const HitRecord& rec, Vec3* attenuation,
               Ray* scattered) const override;

 private:
  // 在 unit sphere 中随机找一个点
  Vec3 random_in_unit_sphere() const {
    Vec3 p;
    do {
      // 随机数是 [0, 1)
      // 但我们需要 (-1, 1)
      p = 2.0 * Vec3(Rand(), Rand(), Rand()) - Vec3(1, 1, 1);
    } while (dot(p, p) >= 1.0);  // 如果随机点不在sphere内,就继续找
    return p;
  }
  Texture *albedo_;
};

bool Lambertian::Scatter(const Ray& ray_in, const HitRecord& rec,
                         Vec3* attenuation, Ray* scattered) const {
  Vec3 target = rec.p + rec.normal + random_in_unit_sphere();
  *scattered = Ray(rec.p, target - rec.p, ray_in.time());
  *attenuation = albedo_->Value(0, 0, rec.p);
  return true;
}

#endif  // __LAMBERTIAN_H__
