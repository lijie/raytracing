// Bounding Volume Hierarchical Tree

#include "hitable.h"
#include "aabb.h"

class BvhNode : public Hitable {
 public:
  BvhNode() {}
  BvhNode(Hitable **list, int size, double t0, double t1);
  bool Hit(const Ray& ray, double t_min, double t_max,
           HitRecord* rec) const override;
  bool BoundingBox(double t0, double t1, AABB* box) const override;
  std::string Name() override { return "BvhNode"; };

 private:
  AABB box_;
  Hitable *l_;
  Hitable *r_;
};

static int bvh_node_compare_x(const void *a, const void *b) {
  AABB box_a, box_b;
  Hitable *ah = *(Hitable **)a;
  Hitable *bh = *(Hitable **)b;

  // printf("bvh_node_compare_x: %p, %p\n", ah, bh);

  bool resa = ah->BoundingBox(0, 0, &box_a);
  bool resb = bh->BoundingBox(0, 0, &box_b);
  assert(resa && resb);

  if (box_a.min().x() - box_b.min().x() < 0)
    return -1;
  else
    return 1;
}

static int bvh_node_compare_y(const void *a, const void *b) {
  AABB box_a, box_b;
  Hitable *ah = *(Hitable **)a;
  Hitable *bh = *(Hitable **)b;

  bool resa = ah->BoundingBox(0, 0, &box_a);
  bool resb = bh->BoundingBox(0, 0, &box_b);
  assert(resa && resb);

  if (box_a.min().y() - box_b.min().y() < 0)
    return -1;
  else
    return 1;
}

static int bvh_node_compare_z(const void *a, const void *b) {
  AABB box_a, box_b;
  Hitable *ah = *(Hitable **)a;
  Hitable *bh = *(Hitable **)b;

  bool resa = ah->BoundingBox(0, 0, &box_a);
  bool resb = bh->BoundingBox(0, 0, &box_b);
  assert(resa && resb);

  if (box_a.min().z() - box_b.min().z() < 0)
    return -1;
  else
    return 1;
}

// 构造 bvh tree
BvhNode::BvhNode(Hitable **list, int size, double t0, double t1) {
  // 随机选一个 axis 作为排序依据
  int axis = int(3 * Rand());

  // printf("BvhNode: %p, %p\n", list[0], list[1]);

  if (axis == 0) {
    // axis-x
    qsort(list, size, sizeof(Hitable *), bvh_node_compare_x);
  } else if (axis == 1) {
    qsort(list, size, sizeof(Hitable *), bvh_node_compare_y);
  } else {
    qsort(list, size, sizeof(Hitable *), bvh_node_compare_z);
  }

  // fast process
  if (size == 1) {
    l_ = r_ = list[0];
  } else if (size == 2) {
    l_ = list[0];
    r_ = list[1];
  } else {
    // 递归构造
    l_ = new BvhNode(list, size / 2, t0, t1);
    r_ = new BvhNode(list + size / 2, size - size / 2, t0, t1);
  }

  // 递归构造box
  AABB left_box, right_box;
  if (!l_->BoundingBox(t0, t1, &left_box) || !r_->BoundingBox(t0, t1, &right_box)) {
    assert(0);
  }
  box_ = AABB::SurroudingBox(left_box, right_box);
}

bool BvhNode::Hit(const Ray& ray, double t_min, double t_max,
                  HitRecord* rec) const {
  if (!box_.Hit(ray, t_min, t_max))
      return false;

  HitRecord l_rec, r_rec;

  bool hit_l = l_ != NULL && l_->Hit(ray, t_min, t_max, &l_rec);
  bool hit_r = r_ != NULL && r_->Hit(ray, t_min, t_max, &r_rec);

  if (hit_l && hit_r) {
    *rec = l_rec.t <= r_rec.t ? l_rec : r_rec;
    return true;
  }
  if (hit_l) {
    *rec = l_rec;
    return true;
  }
  if (hit_r) {
    *rec = r_rec;
    return true;
  }
  return false;
}

bool BvhNode::BoundingBox(double t0, double t1, AABB* box) const {
  *box = box_;
  return true;
}
