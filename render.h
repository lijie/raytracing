#include <thread>
#include <vector>
#include "camera.h"
#include "hitable.h"
#include "material.h"
#include "rand.h"

class Renderer {
 private:
  int nr_of_thread_ = 8;

 public:
  Renderer(int nr_of_thread = 8): nr_of_thread_(nr_of_thread) {}

  // depth 控制反射次数
  virtual Vec3 color_of_ray(const Ray& ray, const Hitable* world, int depth);

  // Render entire world
  void Render(int nx, int ny, const Camera& camera, const Hitable* world,
              std::vector<int>* data);

};

// depth 控制反射次数
Vec3 Renderer::color_of_ray(const Ray& ray, const Hitable* world, int depth) {
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

// Render entire world
void Renderer::Render(int nx, int ny, const Camera& camera,
                      const Hitable* world, std::vector<int>* data) {
  const int thread_count = nr_of_thread_;
  int ns = 100;  // for antialiasing

  data->resize(nx * ny * 3);

  printf("total size:%d\n", data->size());

  auto worker = [](int nx, int ny, int ns, int start_y, int count,
                   const Camera& camera, const Hitable* world, std::vector<int>* data,
                   int offset, Renderer* env) {
    // 遍历像素点, PPM 定义的像素起始点为左上角, 所以从 ny-1 开始
    printf("size: %dx%d, worker range %d, %d, offset: %d\n", nx, ny, start_y,
           start_y - count, offset);
    for (int i = start_y; i > start_y - count; i--) {
      for (int j = 0; j < nx; j++) {
        Vec3 color(0, 0, 0);
        // 抗锯齿, 随机ns次与附近的颜色平均
        for (int s = 0; s < ns; s++) {
          // 计算当前像素点的uv (相对于左下角的偏移)
          double u = double(j + Rand()) / double(nx);
          double v = double(i + Rand()) / double(ny);
          auto ray = camera.GetRay(u, v);
          color += env->color_of_ray(ray, world, 0);
        }
        color /= double(ns);
        (*data)[offset++] = color.x();
        (*data)[offset++] = color.y();
        (*data)[offset++] = color.z();
      }
    }
  };

  int ny_of_per_thread = ny / thread_count;
  int rest = ny - (ny_of_per_thread * thread_count);
  std::vector<std::thread*> vector_of_thread;

  for (int i = 0; i < thread_count; i++) {
    int start_y = ny - i * ny_of_per_thread;
    int count = ny_of_per_thread;
    if (i == thread_count - 1) {
      count += rest;
    }
    int offset = i * ny_of_per_thread * nx * 3;
    std::thread* thread = new std::thread(worker, nx, ny, ns, start_y, count,
                                          camera, world, data, offset, this);
    vector_of_thread.push_back(thread);
    printf("create thread %d\n", i);
  }

  for (int i = 0; i < thread_count; i++) {
    std::thread* thread = vector_of_thread[i];
    thread->join();
  }
}

#if 0
void Scan(int nx, int ny, const Camera& camera, const char* output) {
    std::vector<int> data;
    int ns = 100;  // for antialiasing

    Hitable* world = CreateWorld();

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
          color += color_of_ray(ray, world, 0);
        }
        color /= double(ns);
        data.push_back(color.x());
        data.push_back(color.y());
        data.push_back(color.z());
      }
    }

    write_ppm(output, nx, ny, data);
    DestroyWorld(world);
  }
#endif