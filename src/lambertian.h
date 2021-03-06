#ifndef __LAMBERTIAN_H__
#define __LAMBERTIAN_H__

#include <memory>

#include "material.h"
#include "pdf.h"
#include "rand.h"
#include "texture.h"

class Lambertian : public Material {
 public:
  Lambertian(const Vec3& a) { albedo_ = new ConstantTexture(a); }
  Lambertian(Texture* tex) : albedo_(tex) {}
  bool Scatter(const Ray& ray_in, const HitRecord& rec, Vec3* attenuation,
               Ray* scattered, double* pdf) const override;
  bool Scatter(const Ray& ray_in, const HitRecord& rec,
               ScatterRecord* srec) const override;
  double ScatteringPdf(const Ray& ray, const HitRecord& rec,
                       const Ray& scattered) const override;

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

  // TODO(explain)
  Vec3 random_unit_vector() const {
    auto a = RandDouble(0, 2 * M_PI);
    auto z = RandDouble(-1, 1);
    auto r = sqrt(1 - z * z);
    return Vec3(r * cos(a), r * sin(a), z);
  }
  Texture* albedo_;
};

bool Lambertian::Scatter(const Ray& ray_in, const HitRecord& rec,
                         Vec3* attenuation, Ray* scattered, double* pdf) const {
  // OrthonormalBasis ob;
  // ob.BuildFromW(rec.normal);

  CosinePdf cosine_pdf(rec.normal);

  auto direction =
      cosine_pdf.Generate();  // unit_vector(ob.Local(RandomCosineDirection()));
  // Vec3 target = rec.p + rec.normal + random_unit_vector();
  // *scattered = Ray(rec.p, target - rec.p, ray_in.time());
  *scattered = Ray(rec.p, direction, ray_in.time());

  // BRDF of Lambertian: BRDF = A * (1 / PI)
  // see pbr-book ch8.3 Lambertian Reflection
  *attenuation = albedo_->Value(rec.u, rec.v, rec.p) / M_PI;
  // TODO(explain)
  // *pdf = dot(rec.normal, unit_vector(scattered->direction())) / M_PI;
  // *pdf = dot(ob.w(), scattered->direction()) / M_PI;
  *pdf = cosine_pdf.Value(direction);
  return true;
}

bool Lambertian::Scatter(const Ray& ray_in, const HitRecord& rec, ScatterRecord* srec) const {
  srec->is_specular = false;
  // BRDF of Lambertian: BRDF = A * (1 / PI)
  // see pbr-book ch8.3 Lambertian Reflection
  srec->attenuation = albedo_->Value(rec.u, rec.v, rec.p) / M_PI;
  srec->pdf = std::make_shared<CosinePdf>(rec.normal);
  return true;
}

double Lambertian::ScatteringPdf(const Ray& ray, const HitRecord& rec,
                                 const Ray& scattered) const {
  auto cosine = dot(rec.normal, unit_vector(scattered.direction()));
  return cosine < 0 ? 0 : cosine / M_PI;
}

#endif  // __LAMBERTIAN_H__
