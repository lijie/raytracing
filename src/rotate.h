#ifndef __ROTATE_H__
#define __ROTATE_H__

#include "aabb.h"
#include "hitable.h"
#include "ray.h"

class Material;

class RotateY : public Hitable {
 public:
  RotateY(Hitable* p, double angle);

  bool Hit(const Ray& ray, double t_min, double t_max,
           HitRecord* rec) const override;
  bool BoundingBox(double t0, double t1, AABB* box) const override;

  std::string Name() override { return "RotateY"; }

 private:
  Hitable* p_;
  double radian_;
  double sin_theta_;
  double cos_theta_;
  AABB bbox_;
  bool has_bbox_ = false;
};

static Vec3 RotateVectorY(const Vec3& origin, double sin_theta,
                          double cos_theta) {
  auto new_x = cos_theta * origin.x() + sin_theta * origin.z();
  auto new_z = -sin_theta * origin.x() + cos_theta * origin.z();
  return Vec3(new_x, origin.y(), new_z);
}

RotateY::RotateY(Hitable* p, double angle) {
  radian_ = angle * M_PI / 180;
  p_ = p;

  sin_theta_ = sin(radian_);
  cos_theta_ = cos(radian_);

  AABB box;
  has_bbox_ = p_->BoundingBox(0, 1, &box);

  Vec3 min(INFINITY, INFINITY, INFINITY);
  Vec3 max(-INFINITY, -INFINITY, -INFINITY);

  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      for (int k = 0; k < 2; k++) {
        auto x = i * box.min().x() + (1 - i) * box.max().x();
        auto y = j * box.min().y() + (1 - j) * box.max().y();
        auto z = k * box.min().z() + (1 - k) * box.max().z();

        printf("x: %f, y: %f, z: %f\n", x, y, z);

        // 旋转公式
        auto new_x = cos_theta_ * x + sin_theta_ * z;
        auto new_z = -sin_theta_ * x + cos_theta_ * z;

        Vec3 tmp(new_x, y, new_z);

        for (int c = 0; c < 3; c++) {
          min[c] = fmin(min[c], tmp[c]);
          max[c] = fmax(max[c], tmp[c]);
        }
      }
    }
  }

  bbox_ = AABB(min, max);
}

bool RotateY::Hit(const Ray& ray, double t_min, double t_max,
                  HitRecord* rec) const {
  auto origin = ray.origin();
  auto direction = ray.direction();

  // 这里旋转射线时, 角度是相反的, 根据三角函数公式
  // cos(-A) = cos(A), sin(-A) = -sin(A)
  // 这里cos_theta_, sin_theta_前面的符号做了变换
  origin[0] = cos_theta_ * ray.origin()[0] - sin_theta_ * ray.origin()[2];
  origin[2] = sin_theta_ * ray.origin()[0] + cos_theta_ * ray.origin()[2];

  direction[0] =
      cos_theta_ * ray.direction()[0] - sin_theta_ * ray.direction()[2];
  direction[2] =
      sin_theta_ * ray.direction()[0] + cos_theta_ * ray.direction()[2];

  Ray rotate_r(origin, direction, ray.time());
  // printf("rotate_r: %f->%f\n", origin.x(), direction.x());

  if (!p_->Hit(rotate_r, t_min, t_max, rec)) return false;

  auto p = RotateVectorY(rec->p, sin_theta_, cos_theta_);
  auto normal = RotateVectorY(rec->normal, sin_theta_, cos_theta_);

  rec->p = p;
  rec->normal = normal;
  return true;
}

bool RotateY::BoundingBox(double t0, double t1, AABB* box) const {
  *box = bbox_;
  return has_bbox_;
}

#endif  // __ROTATE_H__
