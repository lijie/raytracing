#ifndef __TEXTURE_H__
#define __TEXTURE_H__

#include "vec3.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

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

class ImageTexture : public Texture {
 public:
  ImageTexture() {}
  ImageTexture(unsigned char *data, int nx, int ny):data_(data), nx_(nx), ny_(ny) {}
  ImageTexture(const char *path);
  Vec3 Value(double u, double v, const Vec3& point) const override {
    // printf("try get uv: %f, %f\n", u, v);
    int i = u * nx_;
    // 图片从左上角开始，但是v从底部开始
    int j = (1 - v) * ny_ - 0.001;
    if (i < 0) i = 0;
    if (j < 0) j = 0;
    if (i >= nx_) i = nx_ - 1;
    if (j >= ny_) j = ny_ - 1;
    double r = int(data_[nx_ * j * 3 + i * 3 + 0]) / 255.0;
    double g = int(data_[nx_ * j * 3 + i * 3 + 1]) / 255.0;
    double b = int(data_[nx_ * j * 3 + i * 3 + 2]) / 255.0;
    return Vec3(r, g, b);
  };

  unsigned char *data_;
  int nx_, ny_;
};

ImageTexture::ImageTexture(const char *path) {
  int comp;
  data_ = stbi_load(path, &nx_, &ny_, &comp, 3);
  if (data_ == NULL) {
    assert(0);
  }
  printf("load image: %s, %dx%d, channel:%d\n", path, nx_, ny_, comp);
}

#endif  // __TEXTURE_H__
