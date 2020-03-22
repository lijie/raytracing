#ifndef __DIELECTRIC_H__
#define __DIELECTRIC_H__

#include "material.h"
#include "rand.h"

class Dielectric : public Material {
 public:
  Dielectric(double ref_idx) : ref_idx_(ref_idx) {}
  bool Scatter(const Ray& ray_in, const HitRecord& rec, Vec3* attenuation,
               Ray* scattered) const override;
  Vec3 reflect(const Vec3& v, const Vec3& n) const;

 private:
  bool refract(const Vec3& v, const Vec3& n, double ni_over_nt,
               Vec3* refracted) const;
  double schlick(double cosine, double ref_idx) const;
  // 参考折射率
  // 空气:1, 水:1.33-1.34, 窗户玻璃: 1.51, 光学玻璃: 1.49-1.92, 钻石: 2.42
  double ref_idx_;
};

// 折射的数学原理参考 FOCG, p325
// tips:
// To add transparent materials to our code, we need a way to determine when
// a ray is going “into” an object. The simplest way to do this is to assume
// that all objects are embedded in air with refractive index very close to 1.0,
// and that surface normals point “out” (toward the air).
bool Dielectric::refract(const Vec3& v, const Vec3& n, double ni_over_nt,
                         Vec3* refracted) const {
  auto unit_v = unit_vector(v);
  auto dt = dot(unit_v, n);
  auto discriminat = 1 - (ni_over_nt * ni_over_nt * (1 - dt * dt));
  if (discriminat <= 0) return false;
  *refracted = ni_over_nt * (unit_v - n * dt) / 1 - n * sqrt(discriminat);
  return true;
}

// see FOCG, p238
Vec3 Dielectric::reflect(const Vec3& v, const Vec3& n) const {
  // 注意这里 dot(v, n) 主要是用来计算 cos(theta), theta 是v与n的夹角
  // dot(v, n) = ||v|| * ||n|| * cos(theta) = cos(theta), 单位向量模都是1
  // 详细请参考教科书的图解
  return v - 2 * dot(v, n) * n;
}

double Dielectric::schlick(double cosine, double ref_idx) const {
  double r0 = (1 - ref_idx) / (1 + ref_idx);
  r0 = r0 * r0;
  return r0 + (1 - r0) * pow(1 - cosine, 5);
}

bool Dielectric::Scatter(const Ray& ray_in, const HitRecord& rec,
                         Vec3* attenuation, Ray* scattered) const {
  Vec3 outward_normal;
  Vec3 reflected = reflect(ray_in.direction(), rec.normal);
  Vec3 refracted;
  double ni_over_nt;
  double reflect_prob;
  double cosine;

  *attenuation = Vec3(1.0, 1.0, 1.0);  // 全反射

  // 这里我们总是假设材质外是空气， 折射率是1， 简化了计算
  if (dot(ray_in.direction(), rec.normal) > 0) {
    // 材质表面由外至内发生折射
    outward_normal = -rec.normal;
    ni_over_nt = ref_idx_;
    // ni_over_nt = 1 / ref_idx_;
    cosine = ref_idx_ * dot(ray_in.direction(), rec.normal) /
             ray_in.direction().length();
  } else {
    // 材质由内至外发生折射
    // 由于外部是空气，折射率是1， 根据 Snell's Law: 入射介质折射率*sin(入射角)
    // = 折射介质折射率*sin(折射角) 两边都除以入射介质折射率:
    outward_normal = rec.normal;
    // ni_over_nt = ref_idx_;
    ni_over_nt = 1 / ref_idx_;
    cosine = -dot(ray_in.direction(), rec.normal) / ray_in.direction().length();
  }

  // 计算折射射线
  if (refract(ray_in.direction(), outward_normal, ni_over_nt, &refracted)) {
    reflect_prob = schlick(cosine, ref_idx_);
    if (Rand() >= reflect_prob) {
      *scattered = Ray(rec.p, refracted);
      return true;
    }
  }
  *scattered = Ray(rec.p, reflected);
  return true;
}

#endif  // #define __DIELECTRIC_H__
