#ifndef __SCENES_RANDOMWORLD_H__
#define __SCENES_RANDOMWORLD_H__

#include "camera.h"
#include "dielectric.h"
#include "lambertian.h"
#include "metal.h"
#include "moving_sphere.h"
#include "scene.h"
#include "sphere.h"
#include "vec3.h"

class TestRandomWorld : public Scene {
 public:
  virtual Hitable* CreateDiffuseSphere(const Vec3& center, double r) {
    auto a = Vec3(Rand() * Rand(), Rand() * Rand(), Rand() * Rand());
    // printf("world_add(sphere_new(vec3(%f, %f, %f), %f,
    // material_new_lambertian(vec3(%f, %f, %f))));\n", center.x(), center.y(),
    // center.z(), r, a.x(), a.y(), a.z());
    return new MovingSphere(center, center + Vec3(0, 0.5 * Rand(), 0), 0, 1, r,
                            new Lambertian(a), "MovingSphere");
  }

  virtual Hitable* CreateBaseSphere(const Vec3 center, double r) {
    return new Sphere(center, r, new Lambertian(Vec3(0.5, 0.5, 0.5)));
  }

  virtual Material* CreateMainSphereMaterial1() { return new Dielectric(1.5); }
  virtual Material* CreateMainSphereMaterial2() {
    return new Lambertian(Vec3(0.4, 0.2, 0.1));
  }
  virtual Material* CreateMainSphereMaterial3() {
    return new Metal(Vec3(0.7, 0.6, 0.5), 0.0);
  }

  void SetupWorld() {
    int n = 500;
    // Hitable** list = new Hitable*[n + 1];
    AddHitable(CreateBaseSphere(Vec3(0, -1000, 0), 1000));
    int i = 1;
    for (int a = -11; a < 11; a++) {
      for (int b = -11; b < 11; b++) {
        float choose_mat = Rand();
        Vec3 center(a + 0.9 * Rand(), 0.2, b + 0.9 * Rand());
        if ((center - Vec3(4, 0.2, 0)).length() > 0.9) {
          if (choose_mat < 0.8) {  // diffuse
#if 0
                new Sphere(center, 0.2,
                           new Lambertian(Vec3(Rand() * Rand(), Rand() * Rand(),
                                               Rand() * Rand())));
#else
            AddHitable(CreateDiffuseSphere(center, 0.2));
#endif
          } else if (choose_mat < 0.95) {  // metal
            auto albedo = Vec3(0.5 * (1 + Rand()), 0.5 * (1 + Rand()),
                               0.5 * (1 + Rand()));
            double fuzz = 0.5 * Rand();
            // printf("world_add(sphere_new(vec3(%f, %f, %f), %f,
            // material_new_mental(vec3(%f, %f, %f), %f)));\n", center.x(),
            // center.y(), center.z(), 0.2, albedo.x(), albedo.y(), albedo.z(),
            // fuzz);
            AddHitable(new Sphere(center, 0.2, new Metal(albedo, fuzz)));
          } else {  // glass
                    // printf("world_add(sphere_new(vec3(%f, %f, %f), %f,
            // material_new_dielectric(%f)));\n", center.x(), center.y(),
            // center.z(), 0.2, 1.5);
            AddHitable(new Sphere(center, 0.2, new Dielectric(1.5)));
          }
        }
      }
    }

    AddHitable(new Sphere(Vec3(0, 1, 0), 1.0, CreateMainSphereMaterial1()));
    AddHitable(new Sphere(Vec3(-4, 1, 0), 1.0, CreateMainSphereMaterial2()));
    AddHitable(new Sphere(Vec3(4, 1, 0), 1.0, CreateMainSphereMaterial3()));

    // return new HitableList(list, i);
    // return new BvhNode(list, i, 0, 1);
  }

  void SetupCamera() {
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
    Camera* camera2 = new Camera(Vec3(13, 2, 3), Vec3(0, 0, 0), Vec3(0, 1, 0),
                                 20, 2, aperture, dist_to_focus, 0, 1);
    SetCamera(camera2);
  }

  void SetupScreen() { SetScreenSize(800, 400); }

  void OnCreate() override {
    SetOutput("test_random_world_2.ppm");
    SetupWorld();
    SetupCamera();
    SetupScreen();
  }
};

#endif  // __SCENES_RANDOMWORLD_H__
