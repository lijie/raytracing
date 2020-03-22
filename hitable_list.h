#ifndef __HITABLE_LIST_H__
#define __HITABLE_LIST_H__

#include "hitable.h"

class HitableList : public Hitable {
 public:
  HitableList() {}
  HitableList(Hitable** list, int size) : list_(list), size_(size) {}
  bool Hit(const Ray& ray, double t_min, double t_max,
           HitRecord* rec) const override;
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

#endif  // __HITABLE_LIST_H__
