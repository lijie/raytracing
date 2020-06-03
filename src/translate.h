#ifndef __TRANSLATE_H__
#define __TRANSLATE_H__

#include "aabb.h"
#include "hitable.h"
#include "ray.h"

class Material;

class Translate : public Hitable {
 public:
  Translate(Hitable* p, const Vec3& offset) : p_(p), offset_(offset) {}

  bool Hit(const Ray& ray, double t_min, double t_max,
           HitRecord* rec) const override;
  bool BoundingBox(double t0, double t1, AABB* box) const override;

  std::string Name() override { return "Translate"; }

 private:
  Hitable* p_;
  Vec3 offset_;
};

bool Translate::Hit(const Ray& ray, double t_min, double t_max,
                    HitRecord* rec) const {
  Ray move_ray(ray.origin() - offset_, ray.direction(), ray.time());
  if (!p_->Hit(move_ray, t_min, t_max, rec)) return false;
  rec->p = rec->p + offset_;
  return true;
}

bool Translate::BoundingBox(double t0, double t1, AABB* box) const {
  if (!p_->BoundingBox(t0, t1, box)) return false;
  *box = AABB(box->min() + offset_, box->max() + offset_);
  return true;
}

#endif  // __TRANSLATE_H__
