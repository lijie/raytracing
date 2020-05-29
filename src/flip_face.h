#ifndef __SRC_FLIP_FACE_H__
#define __SRC_FLIP_FACE_H__

#include "hitable.h"

class FlipFace : public Hitable {
 public:
  FlipFace(Hitable* p) : p_(p) {}

  virtual bool Hit(const Ray& ray, double t_min, double t_max,
                   HitRecord* rec) const {
    if (!p_->Hit(ray, t_min, t_max, rec)) return false;

    rec->normal = -rec->normal;
    return true;
  }

  virtual bool BoundingBox(double t0, double t1, AABB* box) const {
      return p_->BoundingBox(t0, t1, box);
  }

  virtual std::string Name()  {
    return "FlipFace";
  }

 private:
  Hitable* p_;
};

#endif  // __SRC_FLIP_FACE_H__
