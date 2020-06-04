#ifndef __METAL_H__
#define __METAL_H__

#include "material.h"
#include "rand.h"

class Metal : public Material {
 public:
  // fuzz: 对反射后的ray加一个随机摆动
  // 个人觉得fuzz为0时更像镜子的反射，加了fuzz之后更像金属的反射
  Metal(const Vec3& a, double fuzz) : albedo_(a), fuzz_(fuzz) {
    if (fuzz_ > 1) fuzz_ = 1;
  }
  bool Scatter(const Ray& ray_in, const HitRecord& rec, Vec3* attenuation,
               Ray* scattered, double *pdf) const override;

  Vec3 reflect(const Vec3& v, const Vec3& n) const;
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

  Vec3 albedo_;
  double fuzz_;
};

bool Metal::Scatter(const Ray& ray_in, const HitRecord& rec, Vec3* attenuation,
                    Ray* scattered, double *pdf) const {
  auto reflected = reflect(unit_vector(ray_in.direction()), rec.normal);
  *scattered = Ray(rec.p, reflected + fuzz_ * random_in_unit_sphere());
  *attenuation = albedo_;
  // 我的理解是 dot(a,b) 可以看作a投影到b的length
  // 如果反射的射线在法线上的投影<=0, 那说明从投射方向看不到这个反射.
  *pdf = 1.0;
  return (dot(scattered->direction(), rec.normal) > 0);
}

// see FOCG, p238
Vec3 Metal::reflect(const Vec3& v, const Vec3& n) const {
  // 注意这里 dot(v, n) 主要是用来计算 cos(theta), theta 是v与n的夹角
  // dot(v, n) = ||v|| * ||n|| * cos(theta) = cos(theta), 单位向量模都是1
  // 详细请参考教科书的图解
  return v - 2 * dot(v, n) * n;
}

#endif  // __METAL_H__
