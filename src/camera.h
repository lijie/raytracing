#ifndef __CAMERA_H__
#define __CAMERA_H__

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "math.h"
#include "vec3.h"
#include "ray.h"
#include "rand.h"

class Camera {
 public:
  Camera() {
    origin_ = Vec3(0, 0, 0);
    horizontal_ = Vec3(4, 0, 0);
    vertical_ = Vec3(0, 2, 0);
    lower_left_corner_ = Vec3(-2, -1, -1);
  }

  // fov: field of view, 视野, 角度
  // aspect: 宽高比
  Camera(double vfov, double aspect) {
    // Camera 相关的数学原理请参考: FOCG, p144, 7.1.3 The Camera Transformation
    double theta = vfov * M_PI / 180;
    double half_height = tan(theta / 2);
    double half_width = aspect * half_height;

    origin_ = Vec3(0, 0, 0);
    horizontal_ = Vec3(half_width * 2, 0, 0);
    vertical_ = Vec3(0, half_height * 2, 0);
    lower_left_corner_ = Vec3(-half_width, -half_height, -1.0);
  }

  Camera(const Vec3& lookfrom, const Vec3& lookat, Vec3 vup, double vfov,
         double aspect) {
    // Camera 相关的数学原理请参考: FOCG, p144, 7.1.3 The Camera Transformation
    double theta = vfov * M_PI / 180;
    double half_height = tan(theta / 2);
    double half_width = aspect * half_height;

    origin_ = lookfrom;
    // 计算符合右手法则的uvw
    // w: camera 所看方向的反方向
    auto w = unit_vector(lookfrom - lookat);
    // u 必定跟 vup 和 w 垂直
    auto u = unit_vector(cross(vup, w));
    // v 跟w,u垂直
    auto v = cross(w, u);

    lower_left_corner_ = origin_ - half_width * u - half_height * v - w;
    horizontal_ = 2 * half_width * u;
    vertical_ = 2 * half_height * v;
  }

  // Defoucs blur, 散焦模糊??
  // 这部分书上没找到相关内容 :(
  Camera(const Vec3& lookfrom, const Vec3& lookat, Vec3 vup, double vfov,
         double aspect, double aperture, double focus_dist, double t0 = 0,
         double t1 = 0) {
    lens_radius_ = aperture / 2;
    double theta = vfov * M_PI / 180;
    double half_height = tan(theta / 2);
    double half_width = aspect * half_height;

    origin_ = lookfrom;
    // 计算符合右手法则的uvw
    // w: camera 所看方向的反方向
    auto w = unit_vector(lookfrom - lookat);
    // u 必定跟 vup 和 w 垂直
    auto u = unit_vector(cross(vup, w));
    // v 跟w,u垂直
    auto v = cross(w, u);

    lower_left_corner_ = origin_ - half_width * focus_dist * u -
                         half_height * focus_dist * v - focus_dist * w;
    horizontal_ = 2 * half_width * focus_dist * u;
    vertical_ = 2 * half_height * focus_dist * v;

    w_ = w;
    u_ = u;
    v_ = v;

    time0_ = t0;
    time1_ = t1;
  }

  // 指定UV, 返回指向screen的一条射线
  // UV range [0, 1]
  Ray GetRay(double s, double t) const {
    // 生成一个介于 time0 time1 之间的随机时间
    double time = 0;
    if (time0_ != time1_) time = time0_ + Rand() * (time1_ - time0_);
    if (lens_radius_ > 0) {
      Vec3 rd = lens_radius_ * random_in_unit_disk();
      Vec3 offset = u_ * rd.x() + v_ * rd.y();
      return Ray(origin_ + offset,
                 lower_left_corner_ + s * horizontal_ + t * vertical_ -
                     origin_ - offset,
                 time);
    } else {
      return Ray(origin_,
                 lower_left_corner_ + s * horizontal_ + t * vertical_ - origin_,
                 time);
    }
  }

 private:
  Vec3 random_in_unit_disk() const {
    Vec3 p;
    do {
      // 随机数是 [0, 1)
      // 但我们需要 (-1, 1)
      p = 2.0 * Vec3(Rand(), Rand(), 0) - Vec3(1, 1, 0);
    } while (dot(p, p) >= 1.0);  // 如果随机点不在sphere内,就继续找
    return p;
  }
  Vec3 origin_;
  Vec3 lower_left_corner_;
  Vec3 horizontal_;
  Vec3 vertical_;
  Vec3 u_, v_, w_;
  double time0_;
  double time1_;
  double lens_radius_ = 0;
};

#endif  // __CAMERA_H__
