#ifndef __RECT_H__
#define __RECT_H__

#include "aabb.h"
#include "hitable.h"

class Material;
class XYRect : public Hitable {
 public:
  XYRect(double x0, double x1, double y0, double y1, double k, Material *mat)
      : x0_(x0), y0_(y0), x1_(x1), y1_(y1), k_(k), material_(mat) {
    if (x0_ > x1_) std::swap(x0_, x1_);
    if (y0_ > y1_) std::swap(y0_, y1_);
  }
  bool Hit(const Ray &ray, double t_min, double t_max,
           HitRecord *rec) const override;
  bool BoundingBox(double t0, double t1, AABB *box) const override;
  virtual std::string Name() override { return "XYRect"; };

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
  rec->v = (y - y0_) / (y1_ - y0_);
  rec->normal = Vec3(0, 0, 1);
  rec->target = const_cast<XYRect *>(this);
  return true;
}

bool XYRect::BoundingBox(double t0, double t1, AABB *box) const {
  *box = AABB(Vec3(x0_, y0_, k_ - 0.0001), Vec3(x1_, y1_, k_ + 0.0001));
  return true;
}

class XZRect : public Hitable {
 public:
  XZRect(double x0, double x1, double z0, double z1, double k, Material *mat)
      : x0_(x0), x1_(x1), z0_(z0), z1_(z1), k_(k), material_(mat) {
    if (x0_ > x1_) std::swap(x0_, x1_);
    if (z0_ > z1_) std::swap(z0_, z1_);
  }
  bool Hit(const Ray &ray, double t_min, double t_max,
           HitRecord *rec) const override;
  bool BoundingBox(double t0, double t1, AABB *box) const override;
  virtual std::string Name() override { return "XZRect"; };

 private:
  double x0_, x1_, z0_, z1_, k_;
  Material *material_;
};

bool XZRect::Hit(const Ray &ray, double t_min, double t_max,
                 HitRecord *rec) const {
  double t = (k_ - ray.origin().y()) / ray.direction().y();
  if (t < t_min || t > t_max) return false;
  double x = ray.origin().x() + ray.direction().x() * t;
  double z = ray.origin().z() + ray.direction().z() * t;
  if (x < x0_ || x > x1_ || z < z0_ || z > z1_) return false;

  rec->t = t;
  rec->mat = material_;
  rec->p = ray.point_at_parameter(t);
  rec->u = (x - x0_) / (x1_ - x0_);
  rec->v = (z - z0_) / (z1_ - z0_);
  rec->normal = Vec3(0, 1, 0);
  rec->target = const_cast<XZRect *>(this);
  return true;
}

bool XZRect::BoundingBox(double t0, double t1, AABB *box) const {
  *box = AABB(Vec3(x0_, k_ - 0.0001, z0_), Vec3(x1_, k_ + 0.0001, z1_));
  return true;
}


class YZRect : public Hitable {
 public:
  YZRect(double y0, double y1, double z0, double z1, double k, Material *mat)
      : y0_(y0), y1_(y1), z0_(z0), z1_(z1), k_(k), material_(mat) {
    if (y0_ > y1_) std::swap(y0_, y1_);
    if (z0_ > z1_) std::swap(z0_, z1_);
  }
  bool Hit(const Ray &ray, double t_min, double t_max,
           HitRecord *rec) const override;
  bool BoundingBox(double t0, double t1, AABB *box) const override;
  virtual std::string Name() override { return "YZRect"; };

 private:
  double y0_, y1_, z0_, z1_, k_;
  Material *material_;
};

bool YZRect::Hit(const Ray &ray, double t_min, double t_max,
                 HitRecord *rec) const {
  double t = (k_ - ray.origin().x()) / ray.direction().x();
  if (t < t_min || t > t_max) return false;
  double y = ray.origin().y() + ray.direction().y() * t;
  double z = ray.origin().z() + ray.direction().z() * t;
  if (y < y0_ || y > y1_ || z < z0_ || z > z1_) return false;

  rec->t = t;
  rec->mat = material_;
  rec->p = ray.point_at_parameter(t);
  rec->u = (y - y0_) / (y1_ - y0_);
  rec->v = (z - z0_) / (z1_ - z0_);
  rec->normal = Vec3(1, 0, 0);
  rec->target = const_cast<YZRect *>(this);
  return true;
}

bool YZRect::BoundingBox(double t0, double t1, AABB *box) const {
  *box = AABB(Vec3(k_ - 0.0001, y0_, z0_), Vec3(k_ + 0.0001, y1_, z1_));
  return true;
}


#endif  // __RECT_H__