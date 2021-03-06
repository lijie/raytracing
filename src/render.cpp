#include "render.h"

#include <thread>

#include "camera.h"
#include "hitable.h"
#include "material.h"
#include "rand.h"
#include "scene.h"
#include "pdf.h"

// depth 控制反射次数
Vec3 Renderer::color_of_ray(const Ray& ray, const Scene* scene, int depth,
                            bool trace) {
  HitRecord rec;
  ScatterRecord srec;
  Hitable* world = scene->GetRoot();

  if (depth == 0) return Vec3(0, 0, 0);

  if (world->Hit(ray, 0.001, __FLT_MAX__, &rec)) {
    Vec3 emited = rec.mat->Emmited(rec.u, rec.v, rec.p);

    // 如果打到光源, 直接结束
    if (rec.target == scene->GetLightSource())
      return emited;

    if (rec.mat->Scatter(ray, rec, &srec)) {
      auto light_source = scene->GetLightSource();
      double ratio = 0.5;

      if (srec.is_specular) {
        return srec.attenuation * color_of_ray(srec.specular_ray, scene, depth - 1, trace);
      }

      ShapePdf shape_pdf(light_source, rec.p);
      MixturePdf mix_pdf(srec.pdf.get(), &shape_pdf, ratio);

      auto direction = mix_pdf.Generate();
      auto pdf = mix_pdf.Value(direction);
      auto scattered = Ray(rec.p, direction, ray.time());

      auto cosine_theta = dot(rec.normal, unit_vector(scattered.direction()));
      if (cosine_theta < 0)
        cosine_theta = 0;

      // see The Renderring Equation
      // attenuation is value of material BxDF
      return emited + srec.attenuation * cosine_theta * color_of_ray(scattered, scene, depth - 1, trace) / pdf;
    } else {
      return emited;
    }
  } else {
    return scene->BackgroundColor(ray);
  }
}

static void GammaCorrection(int samples_per_pixel, Vec3* color) {
  auto scale = 1.0 / samples_per_pixel;

  auto r = sqrt(scale * color->r());
  auto g = sqrt(scale * color->g());
  auto b = sqrt(scale * color->b());

  *color = Vec3(r, g, b);
}

// Render entire world
void Renderer::Render(int nx, int ny, const Camera& camera, const Scene* scene,
                      std::vector<int>* data) {
  const int thread_count = nr_of_thread_;
  int ns = 50;  // for antialiasing

  data->resize(nx * ny * 3);

  // printf("total size:%d\n", data->size());

  auto worker = [](int nx, int ny, int ns, int start_y, int count,
                   const Camera& camera, const Scene* scene,
                   std::vector<int>* data, int offset, Renderer* env) {
    // 遍历像素点, PPM 定义的像素起始点为左上角, 所以从 ny-1 开始
    // printf("size: %dx%d, worker range %d, %d, offset: %d\n", nx, ny, start_y,
    //        start_y - count, offset);
    for (int i = start_y; i > start_y - count; i--) {
      for (int j = 0; j < nx; j++) {
        Vec3 color(0, 0, 0);
        // 抗锯齿, 随机ns次与附近的颜色平均
        for (int s = 0; s < ns; s++) {
          // 计算当前像素点的uv (相对于左下角的偏移)
          double u = double(j + Rand()) / double(nx);
          double v = double(i + Rand()) / double(ny);
          auto ray = camera.GetRay(u, v);
          if ((i == 500 - 327) && (j == 280)) {
            color += env->color_of_ray(ray, scene, 50, true);
          } else {
            color += env->color_of_ray(ray, scene, 50, false);
          }
        }
        // color /= double(ns);

        GammaCorrection(ns, &color);

        color = RGB(color);
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
                                          camera, scene, data, offset, this);
    vector_of_thread.push_back(thread);
    // printf("create thread %d\n", i);
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