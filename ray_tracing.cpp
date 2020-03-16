//
// Notes:
// 1. FOCG 指 <Fundamentals of Computer Graphics> 4th, Peter Shirley & Steve
// Marschner

#include <assert.h>
#include <math.h>
#include <time.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <vector>

class Vec3 {
 public:
  Vec3() {}
  Vec3(float e1, float e2, float e3) {
    e[0] = e1;
    e[1] = e2;
    e[2] = e3;
  }

  inline const float x() const { return e[0]; }
  inline const float y() const { return e[1]; }
  inline const float z() const { return e[2]; }

  inline const float r() const { return e[0]; }
  inline const float g() const { return e[1]; }
  inline const float b() const { return e[2]; }

  inline const Vec3& operator+() const { return *this; }
  inline Vec3 operator-() const { return Vec3(-e[0], -e[1], -e[2]); }
  inline float operator[](int i) { return e[i]; }
  // inline float& operator[](int i) { return e[i]; }

  inline Vec3 operator+(const Vec3& v2) {
    return Vec3(e[0] + v2.e[0], e[1] + v2.e[1], e[2] + v2.e[2]);
  }
  inline Vec3 operator-(const Vec3& v2) {
    return Vec3(e[0] - v2.e[0], e[1] - v2.e[1], e[2] - v2.e[2]);
  }
  inline Vec3 operator*(const Vec3& v2) {
    return Vec3(e[0] * v2.e[0], e[1] * v2.e[1], e[2] * v2.e[2]);
  }
  inline Vec3 operator/(const Vec3& v2) {
    return Vec3(e[0] / v2.e[0], e[1] / v2.e[1], e[2] / v2.e[2]);
  }
  inline Vec3 operator/(float t) { return Vec3(e[0] / t, e[1] / t, e[2] / t); }
  inline Vec3 operator*(float t) { return Vec3(e[0] * t, e[1] * t, e[2] * t); }

  inline float dot(const Vec3& v2) {
    return x() * v2.x() + y() * v2.y() + z() * v2.z();
  }

  static inline Vec3 cross(const Vec3& v1, const Vec3& v2) {
    return Vec3(v1.y() * v2.z() - v1.z() * v2.y(),
                -(v1.x() * v2.z() - v1.z() * v2.x()),
                v1.x() * v2.y() - v1.y() * v2.x());
  }

  inline Vec3 cross(const Vec3& v2) {
    return Vec3(y() * v2.z() - z() * v2.y(), -(x() * v2.z() - z() * v2.x()),
                x() * v2.y() - y() * v2.x());
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

  inline Vec3& operator*=(float t);
  inline Vec3& operator/=(float t) {
    e[0] /= t;
    e[1] /= t;
    e[2] /= t;
    return *this;
  }

  inline float length() const {
    return sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
  }
  inline float suqared_length() const {
    return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
  }
  inline void make_unit_vector();

  std::string ToString() const {
    std::stringstream ss;
    ss << "[" << e[0] << "," << e[1] << "," << e[2] << "]";
    return ss.str();
  }

  float e[3];
};

inline Vec3 operator+(const Vec3& v1, const Vec3& v2) {
  return Vec3(v1.x() + v2.x(), v1.y() + v2.y(), v1.z() + v2.z());
}
inline Vec3 operator+(const Vec3& v1, float v) {
  return Vec3(v1.x() + v, v1.y() + v, v1.z() + v);
}
inline Vec3 operator+(float v, const Vec3& v1) {
  return Vec3(v1.x() + v, v1.y() + v, v1.z() + v);
}
inline Vec3 operator*(float t, const Vec3& v) {
  return Vec3(t * v.x(), t * v.y(), t * v.z());
}
inline Vec3 operator/(const Vec3& v, float t) {
  return Vec3(v.x() / t, v.y() / t, v.z() / t);
}
inline Vec3 unit_vector(const Vec3& v) { return v / v.length(); }
inline float dot(const Vec3& v1, const Vec3& v2) {
  return v1.x() * v2.x() + v1.y() * v2.y() + v1.z() * v2.z();
}
// [-1, 1] -> [0, 1]
inline Vec3 normalize(const Vec3& v) { return (v + 1) * 0.5; }
// [0, 1] -> [0, 255]
inline Vec3 RGB(const Vec3& v) { return 255.99 * v; }

// 返回随机数 [0, 1)
float Rand() {
  int v = rand();
  return float(v) / float(RAND_MAX);
}

class Ray {
 public:
  Ray() {}
  Ray(const Vec3& a, const Vec3& b) {
    a_ = a;
    b_ = b;
  }

  std::string ToString() const {
    std::stringstream ss;
    ss << a_.ToString();
    ss << " -> ";
    ss << b_.ToString();
    return ss.str();
  }

  Vec3 origin() const { return a_; }
  Vec3 direction() const { return b_; }
  Vec3 point_at_parameter(float t) const { return a_ + t * b_; }
  Vec3 a_;  // 原点
  Vec3 b_;  // 方向
};

class Hitable;
class Material;

struct HitRecord {
  float t;          // paramter of ray
  Vec3 p;           // ray at point
  Vec3 normal;      // normal vector of this point
  Hitable* target;  // hitted target, for debug
  Material* mat;
};

class Hitable {
 public:
  virtual bool Hit(const Ray& ray, float t_min, float t_max,
                   HitRecord* rec) const = 0;
  virtual std::string Name() = 0;
};

class Sphere : public Hitable {
 public:
  Sphere() {}
  Sphere(Vec3 center, float radius, const std::string& name)
      : center_(center), radius_(radius), name_(name) {}

  bool Hit(const Ray& ray, float t_min, float t_max,
           HitRecord* rec) const override;

  std::string Name() override { return name_; }

  Vec3 center_;
  float radius_;
  std::string name_;
};

// 更完整的检测, 简单验证版本参考 hit_sphere()
// 相关数学公式参考 FOCG, p76&p77
bool Sphere::Hit(const Ray& ray, float t_min, float t_max,
                 HitRecord* rec) const {
  Vec3 oc = ray.origin() - center_;
  float a = dot(ray.direction(), ray.direction());
  float b = 2 * dot(ray.direction(), oc);
  float c = dot(oc, oc) - radius_ * radius_;

  // 射线检测 sphere 的方程是一个二元一次方程, discriminat 用来判定方程是否有解
  // see FOCG, p77
  float discriminat = b * b - 4 * a * c;
  // imaginary
  if (discriminat <= 0) return false;

  float t = (-b - sqrt(discriminat)) / (2 * a);
  if (t <= t_min || t >= t_max) {
    t = (-b + sqrt(discriminat)) / (2 * a);
    if (t <= t_min || t >= t_max) {
      return false;
    }
  }

  rec->t = t;
  rec->p = ray.point_at_parameter(t);
  rec->normal = (rec->p - center_) / radius_;
  rec->target = (Hitable*)this;
  // printf("Sphere::Hit: %f\n", t);
  return true;
}

class HitableList : public Hitable {
 public:
  HitableList() {}
  HitableList(Hitable** list, int size) : list_(list), size_(size) {}
  bool Hit(const Ray& ray, float t_min, float t_max,
           HitRecord* rec) const override;

  std::string Name() { return "HitableList"; }

  Hitable** list_;
  int size_;
};

bool HitableList::Hit(const Ray& ray, float t_min, float t_max,
                      HitRecord* rec) const {
  bool hit = false;
  float closest_so_far = t_max;
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

class Camera {
 public:
  Camera() {
    origin_ = Vec3(0, 0, 0);
    horizontal_ = Vec3(4, 0, 0);
    vertical_ = Vec3(0, 2, 0);
    lower_left_corner_ = Vec3(-2, -1, -1);
  }

  // 指定UV, 返回指向screen的一条射线
  // UV range [0, 1]
  Ray GetRay(float u, float v) {
    return Ray(origin_,
               lower_left_corner_ + u * horizontal_ + v * vertical_ - origin_);
  }
  Vec3 origin_;
  Vec3 lower_left_corner_;
  Vec3 horizontal_;
  Vec3 vertical_;
};

// 终于到了材质这一步
class Material {
 public:
  virtual bool Scatter(const Ray& ray_in, const HitRecord& rec,
                       Vec3* attenuation, Ray* scattered) const = 0;
};

class Lambertian : public Material {
 public:
  Lambertian(const Vec3& a) : albedo_(a) {}
  bool Scatter(const Ray& ray_in, const HitRecord& rec, Vec3* attenuation,
               Ray* scattered) const override;

 private:
  // 在 unit sphere 中随机找一个点
  Vec3 random_in_unit_sphere() const {
    Vec3 p;
    do {
      // 随机数是 [0, 1)
      // 但我们需要 (-1, 1)
      p = 2.0 * Vec3(Rand(), Rand(), Rand()) - Vec3(1, 1, 1);
    } while (dot(p, p) >= 1.0);  // 如果随机点不在sphere内,就继续找
    return p;
  }
  Vec3 albedo_;
};

bool Lambertian::Scatter(const Ray& ray_in, const HitRecord& rec,
                         Vec3* attenuation, Ray* scattered) const {
  Vec3 target = rec.p + rec.normal + random_in_unit_sphere();
  *scattered = Ray(rec.p, target - rec.p);
  *attenuation = albedo_;
}

void test_vec3() {
  {
    auto a = Vec3(1, 1, 1);
    auto b = Vec3(2, 2, 2);
    auto c = a + b;
    assert(c.x() == 3 && c.y() == 3 && c.z() == 3);
  }

  {
    auto a = Vec3(3, 4, 5);
    auto b = Vec3(7, 8, 9);
    auto c = a.dot(b);
    assert(c == (21 + 32 + 45));
  }
}

void write_ppm(const char* filename, int nx, int ny,
               const std::vector<int>& data) {
  std::ofstream out;

  out.open(filename);

  out << "P3\n" << nx << " " << ny << "\n255\n";  // PPM header

  for (size_t i = 0; i < data.size(); i++) {
    out << data[i] << "\n";
  }

  out.close();
}

void test_ppm_output() {
  int nx = 200;
  int ny = 100;
  std::vector<int> data;

  for (int i = ny - 1; i >= 0; i--) {
    for (int j = 0; j < nx; j++) {
      float r = float(j) / nx;
      float g = float(i) / ny;
      float b = 0.2;

      int ir = int(255.99 * r);
      int ig = int(255.99 * g);
      int ib = int(255.99 * b);

      // out << ir << "\n" << ig << "\n" << ib << "\n";
      data.push_back(ir);
      data.push_back(ig);
      data.push_back(ib);
    }
  }

  write_ppm("1.ppm", nx, ny, data);
}

// sphere
// formula:
// x*x + y*y + z*z = R*R (在原点)
// (x-cx)*(x-cx) + (y-cy)*(y-cy) + (z-cz)*(z-cz) (在任一点(cx, cy, cz))
// let C = (cx, cy, cz), p = (x, y, z) => dot(p-C, p-C) = R*R
// 把 ray formula P(t) = A + t*B 代入:
// => dot(P(t)-C, P(t)-C) = R*R
// => t*t*dot(B,B)+2*t*dot(B,A-C)+dot(A-C,A-C)-R*R=0 (假定 A 是(0,0,0))
// see FOCG, p76
float hit_sphere(const Vec3& center, float radius, const Ray& ray) {
  Vec3 oc = ray.origin() - center;
  float a = dot(ray.direction(), ray.direction());
  float b = 2 * dot(ray.direction(), oc);
  float c = dot(oc, oc) - radius * radius;

  // 射线检测 sphere 的方程是一个二元一次方程, discriminat 用来判定方程是否有解
  // see FOCG, p77
  float discriminat = b * b - 4 * a * c;

  if (discriminat < 0) {
    // 无解, 说明不相交
    return -1.0;
  } else {
    // 有解的情况下, 取较小的那一个
    return (-b - sqrt(discriminat)) / 2.0 * a;
  }
}

// 右手坐标系
// eye (camera) 在 (0,0,0), screen 4个角 [(-2, 1, -1), (2, 1, -1), (-2, -1, -1),
// (2, -1, -1)] nx, ny 为 screen 的像素尺寸
void linear_blend_blue_to_white(std::vector<int>* out, int nx, int ny) {
  auto color_of_ray = [](const Ray& ray) -> Vec3 {
    auto white = Vec3(1, 1, 1);
    auto blue = Vec3(0.5, 0.7, 1.0);
    auto red = Vec3(1, 0, 0);
    auto unit = unit_vector(ray.direction());

    // 在屏幕中间放一个圆, 检测 ray 是否穿过
    float t = hit_sphere(Vec3(0, 0, -1), 0.5, ray);
    if (t > 0.0) {
      // 有解, 算出法向量, 根据向量减法公式:
      // -(v_center - ray_at_p) -> ray_at_p - v_center
      auto normal = unit_vector(ray.point_at_parameter(t) - Vec3(0, 0, -1));
      // 把法向量转换为rgb显示
      // return 255.99 * ((normal + 1) * 0.5);
      return RGB(normalize(normal));
    }

    // [-1, 1] -> [0, 1]
    float v = (unit.y() + 1.0) * 0.5;
    // blend white & blue
    auto blend_color = (1.0 - v) * white + v * blue;
    // convert to  [0, 255]
    return 255.99 * blend_color;
  };

  // 定义screen的单位尺寸和左下角
  auto lower_left_corner = Vec3(-2, -1, -1);  // 左下角
  auto width = Vec3(4, 0, 0);                 // 宽
  auto height = Vec3(0, 2, 0);                // 高
  // eye所在的位置
  auto origin = Vec3(0, 0, 0);  // 原点

  // 遍历像素点, PPM 定义的像素起始点为左上角, 所以从 ny-1 开始
  for (int i = ny - 1; i >= 0; i--) {
    for (int j = 0; j < nx; j++) {
      // 计算当前像素点的uv (相对于左下角的偏移)
      float u = float(j) / float(nx);
      float v = float(i) / float(ny);
      // 从 eye 发出一条 ray
      // 首先算出 ray 的方向
      auto direction = lower_left_corner + u * width + v * height;
      auto ray = Ray(origin, direction);

      auto color = color_of_ray(ray);

      out->push_back(int(color.r()));
      out->push_back(int(color.g()));
      out->push_back(int(color.b()));
    }
  }
}

void test_linear_blend_blue_to_white() {
  std::vector<int> data;
  int nx = 200;
  int ny = 100;
  linear_blend_blue_to_white(&data, nx, ny);
  write_ppm("test_linear_blend_blue_to_white.ppm", nx, ny, data);
}

void test_two_sphere() {
  std::vector<int> data;
  int nx = 200;
  int ny = 100;
  int ns = 100;  // for antialiasing

  Camera camera;

  // 构造2个球体
  Hitable* list[2];
  list[0] = new Sphere(Vec3(0, 0, -1), 0.5, "sphere_1");
  list[1] = new Sphere(Vec3(0, -100.5, -1), 100, "sphere_2");
  Hitable* world = new HitableList(list, 2);

  auto color_of_ray = [](const Ray& ray, Hitable* world) -> Vec3 {
    HitRecord rec;
    if (world->Hit(ray, 0.0, __FLT_MAX__, &rec)) {
      // 命中的话我们用法向量转RGB作为color
      return RGB(normalize(rec.normal));
    } else {
      // 未命中blend蓝白
      auto white = Vec3(1, 1, 1);
      auto blue = Vec3(0.5, 0.7, 1.0);
      auto unit = unit_vector(ray.direction());
      auto t = (unit.y() + 1.0) * 0.5;
      return RGB((1 - t) * white + t * blue);
    }
  };

  // 遍历像素点, PPM 定义的像素起始点为左上角, 所以从 ny-1 开始
  for (int i = ny - 1; i >= 0; i--) {
    for (int j = 0; j < nx; j++) {
      Vec3 color(0, 0, 0);
      // 抗锯齿, 随机ns次与附近的颜色平均
      for (int s = 0; s < ns; s++) {
        // 计算当前像素点的uv (相对于左下角的偏移)
        float u = float(j + Rand()) / float(nx);
        float v = float(i + Rand()) / float(ny);
        auto ray = camera.GetRay(u, v);
        color += color_of_ray(ray, world);
      }
      color /= float(ns);

      data.push_back(color.x());
      data.push_back(color.y());
      data.push_back(color.z());
    }
  }

  write_ppm("test_two_sphere.ppm", nx, ny, data);

  delete world;
  delete list[0];
  delete list[1];
}

class TestDiffuse {
 private:
  // 在 unit sphere 中随机找一个点
  Vec3 random_in_unit_sphere() {
    Vec3 p;
    do {
      // 随机数是 [0, 1)
      // 但我们需要 (-1, 1)
      p = 2.0 * Vec3(Rand(), Rand(), Rand()) - Vec3(1, 1, 1);
    } while (dot(p, p) >= 1.0);  // 如果随机点不在sphere内,就继续找
    return p;
  }

  Vec3 color_of_ray(const Ray& ray, Hitable* world, int depth = 0) {
    HitRecord rec;
    if (world->Hit(ray, 0.0, __FLT_MAX__, &rec)) {
      // 模拟漫反射, 思路是:
      // 1. 找到 hitpoint, 这里是rec.p, 在p处沿法向量方向构造一个 unit-sphere
      // 2. 在 unit-sphere 内随机寻找一个点 rp,
      // 这里由random_in_unit_sphere()完成
      // 3. 构造一个以 rec.p 为原点, rp方向的射线, 就是反射路径

      // p + normal 是 unit-sphere 的圆心, +rp 就是我们要找的目标, 随机点
      Vec3 rp = random_in_unit_sphere();
      // printf("random point: %s\n", rp.ToString().c_str());
      Vec3 target = rec.p + rec.normal + rp;

      // printf("ray is : %s\n", ray.ToString().c_str());
      // printf("hit: %s, %s, %s\n", rec.target->Name().c_str(),
      // rec.p.ToString().c_str(), rec.normal.ToString().c_str());
      // printf("random target: %s\n", target.ToString().c_str());

      // if (depth > 0) {
      //   printf("hit object %s\n", rec.target->Name().c_str());
      //   printf("target: %s\n", target.ToString().c_str());
      //   exit(0);
      // }
      // 构造新的射线, 方向就是 target, 但要转换到以 rec.p 为原点.
      auto new_ray = Ray(rec.p, target - rec.p);
      // printf("new_ray: %s\n", new_ray.ToString().c_str());
      // 反射光线衰减 50%
      return 0.5 * color_of_ray(new_ray, world, ++depth);
    } else {
      // if (depth > 0) {
      //   printf("color_of_ray, depth > 1, not hit\n");
      // }
      // 未命中blend蓝白
      auto white = Vec3(1, 1, 1);
      auto blue = Vec3(0.5, 0.7, 1.0);
      auto unit = unit_vector(ray.direction());
      auto t = (unit.y() + 1.0) * 0.5;
      return RGB((1 - t) * white + t * blue);
    }
  }

 public:
  void Run() {
    std::cout << "aaaaa" << std::endl;

    std::vector<int> data;
    int nx = 200;
    int ny = 100;
    int ns = 100;  // for antialiasing

    Camera camera;

    // 构造2个球体
    Hitable* list[2];
    list[0] = new Sphere(Vec3(0, 0, -1), 0.5, "sphere_1");
    list[1] = new Sphere(Vec3(0, -100.5, -1), 100, "sphere_2");
    Hitable* world = new HitableList(list, 2);

    // 遍历像素点, PPM 定义的像素起始点为左上角, 所以从 ny-1 开始
    for (int i = ny - 1; i >= 0; i--) {
      for (int j = 0; j < nx; j++) {
        Vec3 color(0, 0, 0);
        // 抗锯齿, 随机ns次与附近的颜色平均
        for (int s = 0; s < ns; s++) {
          // 计算当前像素点的uv (相对于左下角的偏移)
          float u = float(j + Rand()) / float(nx);
          float v = float(i + Rand()) / float(ny);
          auto ray = camera.GetRay(u, v);
          color += color_of_ray(ray, world);
        }
        color /= float(ns);

        data.push_back(color.x());
        data.push_back(color.y());
        data.push_back(color.z());
      }
    }

    std::cout << "aaaaa" << std::endl;
    write_ppm("test_diffuse.ppm", nx, ny, data);

    delete world;
    delete list[0];
    delete list[1];
  }
};

void test_diffuse() {
  TestDiffuse test;
  test.Run();
}

void run_test() {
  // test_vec3();
  // test_ppm_output();
  // test_linear_blend_blue_to_white();
  // test_two_sphere();
  test_diffuse();
}

int main() {
  srand(time(NULL));
  run_test();
}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧:
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5.
//   转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
