#ifndef __BOX_H__
#define __BOX_H__

#include "aabb.h"
#include "hitable.h"
#include "hitable_list.h"
#include "rect.h"
#include "flip_face.h"

class Material;

class Box : public Hitable {
 public:
  Box(Vec3 p0, Vec3 p1, Material* mat);

  bool Hit(const Ray& ray, double t_min, double t_max,
           HitRecord* rec) const override;
  bool BoundingBox(double t0, double t1, AABB* box) const override;

  std::string Name() override { return "Box"; }

 private:
  Vec3 min_;
  Vec3 max_;
  HitableList list;
};

Box::Box(Vec3 p0, Vec3 p1, Material* mat) {
    min_ = p0;
    max_ = p1;

    // 构造6个面
    list.Add(new XYRect(p0.x(), p1.x(), p0.y(), p1.y(), p1.z(), mat));
    list.Add(new FlipFace(new XYRect(p0.x(), p1.x(), p0.y(), p1.y(), p0.z(), mat)));

    // list.Add(new FlipFace(new XZRect(p0.x(), p1.x(), p0.z(), p1.z(), p1.y(), mat)));
    list.Add(new XZRect(p0.x(), p1.x(), p0.z(), p1.z(), p1.y(), mat));
    list.Add(new FlipFace(new XZRect(p0.x(), p1.x(), p0.z(), p1.z(), p0.y(), mat)));

    list.Add(new FlipFace(new YZRect(p0.y(), p1.y(), p0.z(), p1.z(), p0.x(), mat)));
    list.Add(new YZRect(p0.y(), p1.y(), p0.z(), p1.z(), p1.x(), mat));
}

bool Box::Hit(const Ray& ray, double t_min, double t_max,
              HitRecord* rec) const {
  return list.Hit(ray, t_min, t_max, rec);
}

bool Box::BoundingBox(double t0, double t1, AABB* box) const {
  *box = AABB(min_, max_);
  return true;
}

#endif  // __BOX_H__
