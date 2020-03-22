#ifndef __MOVING_SPHERE_H__
#define __MOVING_SPHERE_H__

#include "hitable.h"

class Material;

class MovingSphere : public Hitable {
 public:
  MovingSphere(Vec3 center0, Vec3 center1, double t0, double t1, double radius,
               Material* mat, const std::string& name)
      : center0_(center0),
        center1_(center1),
        time0_(t0),
        time1_(t1),
        radius_(radius),
        mat_(mat),
        name_(name) {}

  bool Hit(const Ray& ray, double t_min, double t_max, HitRecord* rec) const;
  std::string Name() override { return name_; }

  // 给定时间计算出 center
  // 由于 Sphere 是移动的, 获取其圆心就必须指定获取哪一时刻的圆心
  Vec3 center(double time) const;

  Vec3 center0_;
  Vec3 center1_;
  double time0_;
  double time1_;
  double radius_;
  Material* mat_;
  std::string name_;
};

Vec3 MovingSphere::center(double time) const {
  return center0_ +
         ((time - time0_) / (time1_ - time0_)) * (center1_ - center0_);
}

bool MovingSphere::Hit(const Ray& ray, double t_min, double t_max,
                       HitRecord* rec) const {
  Vec3 oc = ray.origin() - center(ray.time());
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
  rec->normal = (rec->p - center(ray.time())) / radius_;
  rec->mat = mat_;
  rec->target = (Hitable*)this;
  // printf("Sphere::Hit: %f\n", t);
  return true;
}

#endif  // __MOVING_SPHERE_H__