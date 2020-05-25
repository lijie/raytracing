#ifndef __RECT_H__
#define __RECT_H__

#include "aabb.h"
#include "hitable.h"

class XYRect : public Hitable {
public:
  virtual bool Hit(const Ray &ray, double t_min, double t_max,
                   HitRecord *rec) const = 0;
  virtual bool BoundingBox(double t0, double t1, AABB *box) const = 0;
  virtual std::string Name() = 0;
};

bool XYRect::Hit(const Ray &ray, double t_min, double t_max,
                 HitRecord *rec) const {

                   
                 }

#endif // __RECT_H__