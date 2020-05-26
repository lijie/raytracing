#ifndef __HITABLE_LIST_H__
#define __HITABLE_LIST_H__

#include "hitable.h"
#include "aabb.h"

class HitableList : public Hitable {
 public:
  HitableList() {}
  HitableList(Hitable** list, int size) : list_(list), size_(size) {}
  bool Hit(const Ray& ray, double t_min, double t_max,
           HitRecord* rec) const override;
  bool BoundingBox(double t0, double t1, AABB* box) const;
  int size() { return size_; }
  Hitable** list() { return list_; }

  std::string Name() { return "HitableList"; }

  Hitable** list_;
  int size_;
};

bool HitableList::Hit(const Ray& ray, double t_min, double t_max,
                      HitRecord* rec) const {
  bool hit = false;
  double closest_so_far = t_max;
  HitRecord tmp;
  for (int i = 0; i < size_; i++) {
    if (list_[i]->Hit(ray, t_min, closest_so_far, &tmp)) {
      hit = true;
      closest_so_far = tmp.t;
      *rec = tmp;
    }
  }
  return hit;
}

// 构造一个BoudingBox可以围住List内的所有对象
bool HitableList::BoundingBox(double t0, double t1, AABB* box) const {
  if (size_ == 0) return false;

  if (!list_[0]->BoundingBox(t0, t1, box)) return false;
  AABB tmp;
  for (int i = 1; i < size_; i++) {
    if (!list_[i]->BoundingBox(t0, t1, &tmp)) return false;
    *box = AABB::SurroudingBox(tmp, *box);
  }
  return true;
}

#endif  // __HITABLE_LIST_H__
