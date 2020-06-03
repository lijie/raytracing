// vec3.h

#ifndef __VEC3_H__
#define __VEC3_H__

#include <string>
#include <sstream>
#include <cmath>

class Vec3 {
 public:
  Vec3() {}
  Vec3(double e1, double e2, double e3) {
    e[0] = e1;
    e[1] = e2;
    e[2] = e3;
  }

  inline const double x() const { return e[0]; }
  inline const double y() const { return e[1]; }
  inline const double z() const { return e[2]; }

  inline const double r() const { return e[0]; }
  inline const double g() const { return e[1]; }
  inline const double b() const { return e[2]; }

  inline const Vec3& operator+() const { return *this; }
  inline Vec3 operator-() const { return Vec3(-e[0], -e[1], -e[2]); }
  inline double& operator[](int i) { return e[i]; }
  inline double operator[](int i) const { return e[i]; }
  // inline double& operator[](int i) { return e[i]; }

  inline Vec3 operator+(const Vec3& v2) {
    return Vec3(e[0] + v2.e[0], e[1] + v2.e[1], e[2] + v2.e[2]);
  }
  // inline Vec3 operator-(const Vec3& v2) {
  //   return Vec3(e[0] - v2.e[0], e[1] - v2.e[1], e[2] - v2.e[2]);
  // }
  inline Vec3 operator*(const Vec3& v2) {
    return Vec3(e[0] * v2.e[0], e[1] * v2.e[1], e[2] * v2.e[2]);
  }
  inline Vec3 operator/(const Vec3& v2) {
    return Vec3(e[0] / v2.e[0], e[1] / v2.e[1], e[2] / v2.e[2]);
  }
  // inline Vec3 operator/(double t) { return Vec3(e[0] / t, e[1] / t, e[2] /
  // t); } inline Vec3 operator*(double t) { return Vec3(e[0] * t, e[1] * t,
  // e[2] * t); }

  inline double dot(const Vec3& v2) {
    return x() * v2.x() + y() * v2.y() + z() * v2.z();
  }

  inline Vec3& operator+=(const Vec3& v2) {
    e[0] += v2.x();
    e[1] += v2.y();
    e[2] += v2.z();
    return *this;
  }
  inline Vec3& operator-=(const Vec3& v2);
  inline Vec3& operator*=(const Vec3& v2);
  inline Vec3& operator/=(const Vec3& v2);

  inline Vec3& operator*=(double t);
  inline Vec3& operator/=(double t) {
    e[0] /= t;
    e[1] /= t;
    e[2] /= t;
    return *this;
  }

  inline double length() const {
    return sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
  }
  inline double suqared_length() const {
    return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
  }
  inline void make_unit_vector();

  std::string ToString() const {
    std::stringstream ss;
    ss << "[" << e[0] << "," << e[1] << "," << e[2] << "]";
    return ss.str();
  }

  double e[3];
};

inline Vec3 operator+(const Vec3& v1, const Vec3& v2) {
  return Vec3(v1.x() + v2.x(), v1.y() + v2.y(), v1.z() + v2.z());
}
inline Vec3 operator+(const Vec3& v1, double v) {
  return Vec3(v1.x() + v, v1.y() + v, v1.z() + v);
}
inline Vec3 operator+(double v, const Vec3& v1) {
  return Vec3(v1.x() + v, v1.y() + v, v1.z() + v);
}
inline Vec3 operator-(const Vec3& v1, const Vec3& v2) {
  return Vec3(v1.x() - v2.x(), v1.y() - v2.y(), v1.z() - v2.z());
}
inline Vec3 operator*(double t, const Vec3& v) {
  return Vec3(t * v.x(), t * v.y(), t * v.z());
}
inline Vec3 operator*(const Vec3& v, double t) {
  return Vec3(t * v.x(), t * v.y(), t * v.z());
}
inline Vec3 operator/(const Vec3& v, double t) {
  return Vec3(v.x() / t, v.y() / t, v.z() / t);
}
inline Vec3 unit_vector(const Vec3& v) { return v / v.length(); }
inline double dot(const Vec3& v1, const Vec3& v2) {
  return v1.x() * v2.x() + v1.y() * v2.y() + v1.z() * v2.z();
}
inline Vec3 cross(const Vec3& v1, const Vec3& v2) {
  // return Vec3(v1.y() * v2.z() - v1.z() * v2.y(),
  //             -(v1.x() * v2.z() - v1.z() * v2.x()),
  //             v1.x() * v2.y() - v1.y() * v2.x());
  return Vec3(v1.y() * v2.z() - v1.z() * v2.y(),
              v1.z() * v2.x() - v1.x() * v2.z(),
              v1.x() * v2.y() - v1.y() * v2.x());
}
// [-1, 1] -> [0, 1]
inline Vec3 normalize(const Vec3& v) { return (v + 1) * 0.5; }
// [0, 1] -> [0, 255]
inline Vec3 RGB(const Vec3& v) { return 255.99 * v; }

#endif // __VEC3_H__
