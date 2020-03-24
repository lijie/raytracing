#ifndef __TEXTURE_H__
#define __TEXTURE_H__

#include "vec3.h"

class Texture {
 public:
  virtual Vec3 Value(double u, double v, const Vec3& point) const = 0;
};

// 返回固定颜色的纹理
class ConstantTexture : public Texture {
 public:
  ConstantTexture(Vec3 color): color_(color) {}
  Vec3 Value(double u, double v, const Vec3& point) const override {
    return color_;
  }
 private:
  Vec3 color_;
};

// 画小格子的纹理
class CheckerTexture : public Texture {
 public:
  CheckerTexture(Texture *t0, Texture *t1, double width): even_(t0), odd_(t1), width_(width) {}
  Vec3 Value(double u, double v, const Vec3& point) const override {
    double sines = sin(width_ * point.x()) * sin(width_ * point.y()) * sin(width_ * point.z());
    if (sines < 0)
      return odd_->Value(u, v, point);
    else
      return even_->Value(u, v, point);
  }

 private:
  Texture *even_;
  Texture *odd_;
  double width_;
};

#endif  // __TEXTURE_H__
