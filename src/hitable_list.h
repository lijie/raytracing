#ifndef __HITABLE_LIST_H__
#define __HITABLE_LIST_H__

#include <vector>
#include "hitable.h"
#include "aabb.h"

class HitableList : public Hitable {
 public:
  HitableList() {}
  HitableList(Hitable** list, int size);
  bool Hit(const Ray& ray, double t_min, double t_max,
           HitRecord* rec) const override;
  bool BoundingBox(double t0, double t1, AABB* box) const;
  void Add(Hitable *h) {
    hitable_vecs_.push_back(h);
  }

  std::string Name() { return "HitableList"; }

  std::vector<Hitable *> hitable_vecs_;
};

HitableList::HitableList(Hitable** list, int size) {
  for (int i = 0; i < size; i++) {
    hitable_vecs_.push_back(list[i]);
  }
}

bool HitableList::Hit(const Ray& ray, double t_min, double t_max,
                      HitRecord* rec) const {
  bool hit = false;
  double closest_so_far = t_max;
  HitRecord tmp;
  for (int i = 0; i < hitable_vecs_.size(); i++) {
    if (hitable_vecs_[i]->Hit(ray, t_min, closest_so_far, &tmp)) {
      hit = true;
      closest_so_far = tmp.t;
      *rec = tmp;
    }
  }
  return hit;
}

// 构造一个BoudingBox可以围住List内的所有对象
bool HitableList::BoundingBox(double t0, double t1, AABB* box) const {
  if (hitable_vecs_.size() == 0) return false;

  if (!hitable_vecs_[0]->BoundingBox(t0, t1, box)) return false;
  AABB tmp;
  for (int i = 1; i < hitable_vecs_.size(); i++) {
    if (!hitable_vecs_[i]->BoundingBox(t0, t1, &tmp)) return false;
    *box = AABB::SurroudingBox(tmp, *box);
  }
  return true;
}

#endif  // __HITABLE_LIST_H__
