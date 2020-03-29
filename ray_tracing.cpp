//
// Notes:
// 1. FOCG 指 <Fundamentals of Computer Graphics> 4th, Peter Shirley & Steve
// Marschner

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <assert.h>
#include <math.h>
#include <time.h>

#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <thread>
#include <vector>

#include "bvh.h"
#include "camera.h"
#include "dielectric.h"
#include "hitable.h"
#include "hitable_list.h"
#include "lambertian.h"
#include "material.h"
#include "metal.h"
#include "moving_sphere.h"
#include "rand.h"
#include "ray.h"
#include "render.h"
#include "sphere.h"
#include "vec3.h"

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
      double r = double(j) / nx;
      double g = double(i) / ny;
      double b = 0.2;

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
double hit_sphere(const Vec3& center, double radius, const Ray& ray) {
  Vec3 oc = ray.origin() - center;
  double a = dot(ray.direction(), ray.direction());
  double b = 2 * dot(ray.direction(), oc);
  double c = dot(oc, oc) - radius * radius;

  // 射线检测 sphere 的方程是一个二元一次方程, discriminat 用来判定方程是否有解
  // see FOCG, p77
  double discriminat = b * b - 4 * a * c;

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
    double t = hit_sphere(Vec3(0, 0, -1), 0.5, ray);
    if (t > 0.0) {
      // 有解, 算出法向量, 根据向量减法公式:
      // -(v_center - ray_at_p) -> ray_at_p - v_center
      auto normal = unit_vector(ray.point_at_parameter(t) - Vec3(0, 0, -1));
      // 把法向量转换为rgb显示
      // return 255.99 * ((normal + 1) * 0.5);
      return RGB(normalize(normal));
    }

    // [-1, 1] -> [0, 1]
    double v = (unit.y() + 1.0) * 0.5;
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
      double u = double(j) / double(nx);
      double v = double(i) / double(ny);
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
  list[0] = new Sphere(Vec3(0, 0, -1), 0.5, new Lambertian(Vec3(0.5, 0.5, 0.5)),
                       "sphere_1");
  list[1] = new Sphere(Vec3(0, -100.5, -1), 100,
                       new Lambertian(Vec3(0.5, 0.5, 0.5)), "sphere_2");
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
        double u = double(j + Rand()) / double(nx);
        double v = double(i + Rand()) / double(ny);
        auto ray = camera.GetRay(u, v);
        color += color_of_ray(ray, world);
      }
      color /= double(ns);

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
    std::vector<int> data;
    int nx = 200;
    int ny = 100;
    int ns = 100;  // for antialiasing

    Camera camera;

    // 构造2个球体
    Hitable* list[2];
    list[0] = new Sphere(Vec3(0, 0, -1), 0.5,
                         new Lambertian(Vec3(0.5, 0.5, 0.5)), "sphere_1");
    list[1] = new Sphere(Vec3(0, -100.5, -1), 100,
                         new Lambertian(Vec3(0.5, 0.5, 0.5)), "sphere_2");
    Hitable* world = new HitableList(list, 2);

    // 遍历像素点, PPM 定义的像素起始点为左上角, 所以从 ny-1 开始
    for (int i = ny - 1; i >= 0; i--) {
      for (int j = 0; j < nx; j++) {
        Vec3 color(0, 0, 0);
        // 抗锯齿, 随机ns次与附近的颜色平均
        for (int s = 0; s < ns; s++) {
          // 计算当前像素点的uv (相对于左下角的偏移)
          double u = double(j + Rand()) / double(nx);
          double v = double(i + Rand()) / double(ny);
          auto ray = camera.GetRay(u, v);
          color += color_of_ray(ray, world);
        }
        color /= double(ns);

        data.push_back(color.x());
        data.push_back(color.y());
        data.push_back(color.z());
      }
    }

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

class TestMetal {
 private:
  Camera* camera_ = nullptr;
  const char* output_ = NULL;

  // depth 控制反射次数
  Vec3 color_of_ray(const Ray& ray, Hitable* world, int depth) {
    HitRecord rec;

    if (world->Hit(ray, 0.001, __FLT_MAX__, &rec)) {
      Ray scattered;
      Vec3 attenuation;

      if (depth < 50 && rec.mat->Scatter(ray, rec, &attenuation, &scattered)) {
        return attenuation * color_of_ray(scattered, world, depth + 1);
      } else {
        return Vec3(0, 0, 0);  // black
      }
    } else {
      // 未命中blend蓝白
      auto white = Vec3(1, 1, 1);
      auto blue = Vec3(0.5, 0.7, 1.0);
      auto unit = unit_vector(ray.direction());
      auto t = (unit.y() + 1.0) * 0.5;
      return RGB((1 - t) * white + t * blue);
    }
  }

 public:
  TestMetal(const char* output) {
    output_ = output;
    camera_ = new Camera();
  }
  TestMetal(const char* output, Camera* camera) : camera_(camera) {
    output_ = output;
  }
  void Run() {
    std::vector<int> data;
    int nx = 200;
    int ny = 100;
    int ns = 100;  // for antialiasing

    Hitable* list[5];
    list[0] = new Sphere(Vec3(0, 0, -1), 0.5,
                         new Lambertian(Vec3(0.8, 0.3, 0.3)), "sphere_1");
    list[1] = new Sphere(Vec3(0, -100.5, -1), 100,
                         new Lambertian(Vec3(0.8, 0.8, 0.0)), "sphere_2");
    list[2] = new Sphere(Vec3(1, 0, -1), 0.5,
                         new Metal(Vec3(0.8, 0.6, 0.2), 0.3), "sphere_3");
    // list[3] = new Sphere(Vec3(-1, 0, -1), 0.5, new Metal(Vec3(0.8, 0.8,
    // 0.8), 1.0), "sphere_4");
    list[3] = new Sphere(Vec3(-1, 0, -1), 0.5, new Dielectric(1.5), "sphere_4");
    list[4] =
        new Sphere(Vec3(-1, 0, -1), -0.45, new Dielectric(1.5), "sphere_5");

    Hitable* world = new HitableList(list, 5);
    BvhNode* bvh_world = new BvhNode(list, 5, 0, 0);

    // 遍历像素点, PPM 定义的像素起始点为左上角, 所以从 ny-1 开始
    for (int i = ny - 1; i >= 0; i--) {
      for (int j = 0; j < nx; j++) {
        Vec3 color(0, 0, 0);
        // 抗锯齿, 随机ns次与附近的颜色平均
        for (int s = 0; s < ns; s++) {
          // 计算当前像素点的uv (相对于左下角的偏移)
          double u = double(j + Rand()) / double(nx);
          double v = double(i + Rand()) / double(ny);
          auto ray = camera_->GetRay(u, v);
          color += color_of_ray(ray, bvh_world, 0);
        }
        color /= double(ns);

        data.push_back(color.x());
        data.push_back(color.y());
        data.push_back(color.z());
      }
    }

    write_ppm(output_, nx, ny, data);

    delete world;
    delete list[0];
    delete list[1];
    delete list[2];
    delete list[3];
    delete list[4];
  }
};

void test_metal() {
  TestMetal test("test_metal.ppm");
  test.Run();
}

class BaseTest {
 public:
  // depth 控制反射次数
  virtual Vec3 color_of_ray(const Ray& ray, Hitable* world, int depth) {
    HitRecord rec;

    if (world->Hit(ray, 0.001, __FLT_MAX__, &rec)) {
      Ray scattered;
      Vec3 attenuation;

      if (depth < 50 && rec.mat->Scatter(ray, rec, &attenuation, &scattered)) {
        return attenuation * color_of_ray(scattered, world, depth + 1);
      } else {
        return Vec3(0, 0, 0);  // black
      }
    } else {
      // 未命中blend蓝白
      auto white = Vec3(1, 1, 1);
      auto blue = Vec3(0.5, 0.7, 1.0);
      auto unit = unit_vector(ray.direction());
      auto t = (unit.y() + 1.0) * 0.5;
      return RGB((1 - t) * white + t * blue);
    }
  }

  virtual Hitable* CreateWorld() = 0;
  virtual void DestroyWorld(Hitable* world) = 0;

  virtual void Render(int nx, int ny, const Camera& camera,
                      const char* output) {
    Renderer renderer;
    Hitable* world = CreateWorld();
    std::vector<int> data;
    renderer.Render(nx, ny, camera, world, &data);
    write_ppm(output, nx, ny, data);
    DestroyWorld(world);
  }

  virtual void Run(const char* output) = 0;
};

class TestCamera : public BaseTest {
 public:
  Hitable* CreateWorld() override {
    Hitable** p = new (Hitable* [2]);

    double r = cos(M_PI / 4);
    p[0] = new Sphere(Vec3(-r, 0, -1), r, new Lambertian(Vec3(0, 0, 1)),
                      "sphere_1");
    p[1] = new Sphere(Vec3(r, 0, -1), r, new Lambertian(Vec3(1, 0, 0)),
                      "sphere_2");

    // return new HitableList(p, 2);
    return new BvhNode(p, 2, 0, 0);
  }

  void DestroyWorld(Hitable* world) override {
    // HitableList* list = dynamic_cast<HitableList*>(world);
    // Hitable** plist = list->list();
    // for (int i = 0; i < list->size(); i++) {
    //   delete plist[i];
    // }
    // delete list;
  }

  void Run(const char* output) {
    Camera camera(90, 200 / 200);
    Render(200, 100, camera, "test_camera.ppm");
  }
};

void test_camera() {
  TestCamera test;
  test.Run("");

  auto* camera1 =
      new Camera(Vec3(-2, 2, 1), Vec3(0, 0, -1), Vec3(0, 1, 0), 90, 2);
  TestMetal test1("test_camera2.ppm", camera1);
  test1.Run();

  auto* camera2 =
      new Camera(Vec3(-2, 2, 1), Vec3(0, 0, -1), Vec3(0, 1, 0), 30, 2);
  TestMetal test2("test_camera3.ppm", camera2);
  test2.Run();

  double dist_to_focus = (Vec3(-2, 2, 1) - Vec3(0, 0, -1)).length();
  auto* camera3 = new Camera(Vec3(-2, 2, 1), Vec3(0, 0, -1), Vec3(0, 1, 0), 30,
                             2, 2.0, dist_to_focus);
  TestMetal test3("test_camera4.ppm", camera3);
  test3.Run();
}

class TestRandomWorld : public BaseTest {
 public:
  virtual Hitable* CreateDiffuseSphere(const Vec3& center, double r) {
    return new MovingSphere(
        center, center + Vec3(0, 0.5 * Rand(), 0), 0, 1, r,
        new Lambertian(Vec3(Rand() * Rand(), Rand() * Rand(), Rand() * Rand())),
        "MovingSphere");
  }

  virtual Hitable* CreateBaseSphere(const Vec3 center, double r) {
    return new Sphere(center, r, new Lambertian(Vec3(0.5, 0.5, 0.5)));
  }

  virtual Material* CreateMainSphereMaterial1() {
    return new Dielectric(1.5);
  }
  virtual Material* CreateMainSphereMaterial2() {
    return new Lambertian(Vec3(0.4, 0.2, 0.1));
  }
  virtual Material* CreateMainSphereMaterial3() {
    return new Metal(Vec3(0.7, 0.6, 0.5), 0.0);
  }

  Hitable* CreateWorld() override {
    int n = 500;
    Hitable** list = new Hitable*[n + 1];
    list[0] = CreateBaseSphere(Vec3(0, -1000, 0), 1000);
    int i = 1;
    for (int a = -11; a < 11; a++) {
      for (int b = -11; b < 11; b++) {
        float choose_mat = Rand();
        Vec3 center(a + 0.9 * Rand(), 0.2, b + 0.9 * Rand());
        if ((center - Vec3(4, 0.2, 0)).length() > 0.9) {
          if (choose_mat < 0.8) {  // diffuse
            list[i++] =
#if 0
                new Sphere(center, 0.2,
                           new Lambertian(Vec3(Rand() * Rand(), Rand() * Rand(),
                                               Rand() * Rand())));
#else
                CreateDiffuseSphere(center, 0.2);
#endif
          } else if (choose_mat < 0.95) {  // metal
            list[i++] = new Sphere(
                center, 0.2,
                new Metal(Vec3(0.5 * (1 + Rand()), 0.5 * (1 + Rand()),
                               0.5 * (1 + Rand())),
                          0.5 * Rand()));
          } else {  // glass
            list[i++] = new Sphere(center, 0.2, new Dielectric(1.5));
          }
        }
      }
    }

    list[i++] = new Sphere(Vec3(0, 1, 0), 1.0, CreateMainSphereMaterial1());
    list[i++] =
        new Sphere(Vec3(-4, 1, 0), 1.0, CreateMainSphereMaterial2());
    list[i++] =
        new Sphere(Vec3(4, 1, 0), 1.0, CreateMainSphereMaterial3());

    // return new HitableList(list, i);
    return new BvhNode(list, i, 0, 1);
  }

  void DestroyWorld(Hitable* world) override {}

  void Run(const char* output) {
    Vec3 lookfrom = Vec3(-2, 2, 1);
    Vec3 lookat = Vec3(0, 0, -1);
    double dist_to_focus = (lookfrom - lookat).length();
    double aperture = 0;
#if 0
    Camera camera1(Vec3(-2, 2, 1), Vec3(0, 0, -1), Vec3(0, 1, 0), 90, 2, aperture, dist_to_focus);
    Render(400, 200, camera1, "test_random_world_1.ppm");
#endif
    lookfrom = Vec3(13, 2, 3);
    lookat = Vec3(0, 0, 0);
    dist_to_focus = 10;
    (lookfrom - lookat).length();
    Camera camera2(Vec3(13, 2, 3), Vec3(0, 0, 0), Vec3(0, 1, 0), 20, 2,
                   aperture, dist_to_focus, 0, 1);
    // Render(800, 400, camera2, "test_random_world_2.ppm");
    Render(800, 400, camera2, output);
  }
};

void test_random_world() {
  TestRandomWorld test;
  test.Run("test_random_world_2.ppm");
}

// 改造一下之前的函数测试bvh
void test_two_bvh_1() {
  std::vector<int> data;
  int nx = 200;
  int ny = 100;
  int ns = 100;  // for antialiasing

  Camera camera;

  // 构造2个球体
  Hitable* list[2];
  list[0] = new Sphere(Vec3(0, 0, -1), 0.5, new Lambertian(Vec3(0.5, 0.5, 0.5)),
                       "sphere_1");
  list[1] = new Sphere(Vec3(0, -100.5, -1), 100,
                       new Lambertian(Vec3(0.5, 0.5, 0.5)), "sphere_2");
  Hitable* world = new HitableList(list, 2);
  BvhNode* bvh_world = new BvhNode(list, 2, 0, 0);

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
        double u = double(j + Rand()) / double(nx);
        double v = double(i + Rand()) / double(ny);
        auto ray = camera.GetRay(u, v);
        color += color_of_ray(ray, bvh_world);
      }
      color /= double(ns);

      data.push_back(color.x());
      data.push_back(color.y());
      data.push_back(color.z());
    }
  }

  write_ppm("test_two_sphere_bvh.ppm", nx, ny, data);

  delete world;
  delete list[0];
  delete list[1];
}

class TestRandomWorldWithTexture : public TestRandomWorld {
 public:
  Hitable* CreateBaseSphere(const Vec3 center, double r) override {
    Texture* t0 = new ConstantTexture(Vec3(0.8, 0.2, 0.2));
    Texture* t1 = new ConstantTexture(Vec3(0.9, 0.9, 0.9));
    Texture* tex = new CheckerTexture(t0, t1, 10);
    return new Sphere(center, r, new Lambertian(tex));
  }
};

void test_texture_1() {
  TestRandomWorldWithTexture test;
  test.Run("test_texture_1.ppm");
}

class TestTextureEarth : public BaseTest {
 public:
  Hitable* CreateWorld() override;
  void DestroyWorld(Hitable* world) override;
  void Run(const char* output);
};

Hitable* TestTextureEarth::CreateWorld() {
  auto texture = new ImageTexture("earthmap.jpg");
  auto earth_surface = new Lambertian(texture);
  return new Sphere(Vec3(0, 0, 0), 2, earth_surface, "sphere_1");
}

void TestTextureEarth::DestroyWorld(Hitable* world) {}

void TestTextureEarth::Run(const char* output) {
  Vec3 lookfrom = Vec3(-2, 2, 1);
  Vec3 lookat = Vec3(0, 0, -1);
  double dist_to_focus = (lookfrom - lookat).length();
  double aperture = 0;
  lookfrom = Vec3(13, 2, 3);
  lookat = Vec3(0, 0, 0);
  dist_to_focus = 10;
  (lookfrom - lookat).length();
  Camera camera2(Vec3(13, 2, 3), Vec3(0, 0, 0), Vec3(0, 1, 0), 20, 2, aperture,
                 dist_to_focus, 0, 1);
  Render(400, 200, camera2, output);
}

void test_texture_earth() {
  TestTextureEarth test;
  test.Run("test_texture_earth.ppm");
}

class TestRandomWorldWithEarthTexture : public TestRandomWorld {
  public:
    Material *CreateMainSphereMaterial3() {
      auto texture = new ImageTexture("earthmap.jpg");
      return new Lambertian(texture);
    }
};

void test_random_world_earth_texture() {
  TestRandomWorldWithEarthTexture test;
  test.Run("test_random_world_earth_texture.ppm");
}

void run_test() {
  // test_vec3();
  // test_ppm_output();
  // test_linear_blend_blue_to_white();
  // test_two_sphere();
  // test_diffuse();
  // test_metal();
  // test_camera();
  // test_random_world();
  // test_two_bvh_1();
  // test_texture_1();
  // test_texture_earth();
  test_random_world_earth_texture();
}

int main() {
  srand(time(NULL));
  std::cout << time(NULL) << std::endl;
  run_test();
  std::cout << time(NULL) << std::endl;
}
