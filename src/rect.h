#ifndef __RECT_H__
#define __RECT_H__

#include "aabb.h"
#include "hitable.h"

class Material;
class XYRect : public Hitable {
 public:
  XYRect(double x0, double x1, double y0, double y1, double k, Material *mat)
      : x0_(x0), y0_(y0), x1_(x1), y1_(y1), k_(k), material_(mat) {}
  bool Hit(const Ray &ray, double t_min, double t_max,
           HitRecord *rec) const override;
  bool BoundingBox(double t0, double t1, AABB *box) const override;
  virtual std::string Name() override;

 private:
  double x0_, x1_, y0_, y1_, k_;
  Material *material_;
};

bool XYRect::Hit(const Ray &ray, double t_min, double t_max,
                 HitRecord *rec) const {
  double t = (k_ - ray.origin().z()) / ray.direction().z();
  if (t < t_min || t > t_max) return false;
  double x = ray.origin().x() + ray.direction().x() * t;
  double y = ray.origin().y() + ray.direction().y() * t;
  if (x < x0_ || x > x1_ || y < y0_ || y > y1_) return false;

  rec->t = t;
  rec->mat = material_;
  rec->p = ray.point_at_parameter(t);
  rec->u = (x - x0_) / (x1_ - x0_);
  rec->v = (y - x0_) / (y1_ - y0_);
  rec->normal = Vec3(0, 0, 1);
  // rec->target = this;
  return true;
}

#endif  // __RECT_H__