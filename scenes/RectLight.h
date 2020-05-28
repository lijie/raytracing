#ifndef __SCENE_RECT_LIGHT_H__
#define __SCENE_RECT_LIGHT_H__

#include "camera.h"
#include "diffuse_light.h"
#include "lambertian.h"
#include "rect.h"
#include "scene.h"
#include "sphere.h"
#include "texture.h"

class RectLightScene : public Scene {
 public:
  void SetupCamera() {
    Vec3 lookfrom = Vec3(-2, 2, 1);
    Vec3 lookat = Vec3(0, 0, -1);
    double dist_to_focus = (lookfrom - lookat).length();
    double aperture = 0;
#if 0
    Camera camera1(Vec3(-2, 2, 1), Vec3(0, 0, -1), Vec3(0, 1, 0), 90, 2, aperture, dist_to_focus);
    Render(400, 200, camera1, "RectLightScene.ppm");
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

  void OnCreate() {
    SetupScreen();
    SetupCamera();
    SetOutput("RectLightScene.ppm");
    AddHitable(new Sphere(Vec3(0, -1000, 0), 1000,
                          new Lambertian(Vec3(0.5, 0.5, 0.5))));
    AddHitable(
        new Sphere(Vec3(0, 2, 0), 2, new Lambertian(Vec3(0.5, 0.5, 0.5))));

    auto mat = new DiffuseLight(new ConstantTexture(Vec3(4, 4, 4)));
    AddHitable(new Sphere(Vec3(0, 7, 0), 2, mat));
    AddHitable(new XYRect(3, 5, 1, 3, -2, mat));
  }
};

#endif  // __SCENE_RECT_LIGHT_H__
