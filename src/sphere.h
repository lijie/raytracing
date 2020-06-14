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

  Vec3 SurfaceRandomPoint(const Vec3& origin) override;
  double SurfaceRandomPdf(const Vec3& origin, const Vec3& direction) override;

  // 球形坐标系
  // u = phi / 2pi
  // v = theta / pi
  // phi = atan2(y, x)
  // theta = asin(z)
  static void GetUV(const Vec3& point, double* u, double* v) {
    double phi = atan2(point.z(), point.x());  // atan2 返回 [-pi, pi]
    double theta = asin(point.y());            // asin 返回 [-pi/2, pi/2]
    *u = 1 - (phi + M_PI) / (2 * M_PI);
    *v = (theta + M_PI / 2) / M_PI;
  }

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
  GetUV((rec->p - center_) / radius_, &rec->u, &rec->v);
  // printf("Sphere::Hit: %f\n", t);
  return true;
}

bool Sphere::BoundingBox(double t0, double t1, AABB* box) const {
  *box = AABB(center_ - Vec3(radius_, radius_, radius_),
              center_ + Vec3(radius_, radius_, radius_));
  return true;
}

Vec3 Sphere::SurfaceRandomPoint(const Vec3& origin) {
  Vec3 direction = center_ - origin;
  auto distance_squared = direction.suqared_length();
  OrthonormalBasis ob;
  ob.BuildFromW(direction);
  return ob.Local(RandomToSphere(radius_, distance_squared));
}

double Sphere::SurfaceRandomPdf(const Vec3& origin, const Vec3& direction) {
  HitRecord rec;
  if (!this->Hit(Ray(origin, direction), 0.0001, __FLT_MAX__, &rec)) return 0;

  auto cos_theta_max =
      sqrt(1 - radius_ * radius_ / (center_ - origin).suqared_length());
  auto solid_angle = 2 * M_PI * (1 - cos_theta_max);

  return 1 / solid_angle;
}

#endif  // __SPHERE_H__
