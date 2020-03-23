#ifndef __SPHERE_H__
#define __SPHERE_H__

#include "aabb.h"
#include "hitable.h"

class Material;

class Sphere : public Hitable {
 public:
  Sphere() {}
  Sphere(Vec3 center, double radius, Material* mat, const std::string& name)
      : center_(center), radius_(radius), mat_(mat), name_(name) {}
  Sphere(Vec3 center, double radius, Material* mat)
      : center_(center), radius_(radius), mat_(mat) {
    name_ = "unknowm sphere";
  }

  bool Hit(const Ray& ray, double t_min, double t_max,
           HitRecord* rec) const override;
  bool BoundingBox(double t0, double t1, AABB* box) const override;

  std::string Name() override { return name_; }

  Vec3 center_;
  double radius_;
  Material* mat_;
  std::string name_;
};

// 更完整的检测, 简单验证版本参考 hit_sphere()
// 相关数学公式参考 FOCG, p76&p77
bool Sphere::Hit(const Ray& ray, double t_min, double t_max,
                 HitRecord* rec) const {
  Vec3 oc = ray.origin() - center_;
  double a = dot(ray.direction(), ray.direction());
  double b = 2 * dot(ray.direction(), oc);
  double c = dot(oc, oc) - radius_ * radius_;

  // 射线检测 sphere 的方程是一个二元一次方程, discriminat 用来判定方程是否有解
  // see FOCG, p77
  double discriminat = b * b - 4 * a * c;
  // imaginary
  if (discriminat <= 0) return false;

  double t = (-b - sqrt(discriminat)) / (2 * a);
  if (t <= t_min || t >= t_max) {
    t = (-b + sqrt(discriminat)) / (2 * a);
    if (t <= t_min || t >= t_max) {
      return false;
    }
  }

  rec->t = t;
  rec->p = ray.point_at_parameter(t);
  rec->normal = (rec->p - center_) / radius_;
  rec->mat = mat_;
  rec->target = (Hitable*)this;
  // printf("Sphere::Hit: %f\n", t);
  return true;
}

bool Sphere::BoundingBox(double t0, double t1, AABB* box) const {
  *box = AABB(center_ - Vec3(radius_, radius_, radius_),
              center_ + Vec3(radius_, radius_, radius_));
  return true;
}

#endif  // __SPHERE_H__