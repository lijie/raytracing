#ifndef __SCENE_CORNELL_BOX_H__
#define __SCENE_CORNELL_BOX_H__

#include "box.h"
#include "camera.h"
#include "diffuse_light.h"
#include "flip_face.h"
#include "lambertian.h"
#include "rect.h"
#include "scene.h"
#include "sphere.h"
#include "texture.h"

class CornellBox : public Scene {
 public:
  const int width = 1000;
  const int height = 1000;
  void SetupCamera() {
    Vec3 lookfrom = Vec3(278, 278, -800);
    Vec3 lookat = Vec3(278, 278, 0);
    double dist_to_focus = 10;
    Camera* camera2 = new Camera(lookfrom, lookat, Vec3(0, 1, 0), 40,
                                 width / height, 0, dist_to_focus, 0, 1);
    SetCamera(camera2);
  }

  void SetupScreen() { SetScreenSize(width, height); }

  void OnCreate() {
    SetupScreen();
    SetupCamera();
    SetOutput("CornellBox.ppm");

    auto red_mat = new Lambertian(new ConstantTexture(Vec3(0.65, 0.05, 0.05)));
    auto white_mat =
        new Lambertian(new ConstantTexture(Vec3(0.73, 0.73, 0.73)));
    auto green_mat =
        new Lambertian(new ConstantTexture(Vec3(0.12, 0.45, 0.15)));
    auto light_mat = new DiffuseLight(new ConstantTexture(Vec3(15, 15, 15)));

    AddHitable(new FlipFace(new YZRect(0, 555, 0, 555, 555, green_mat)));
    AddHitable(new YZRect(0, 555, 0, 555, 0, red_mat));
    AddHitable(new XZRect(213, 343, 227, 332, 554, light_mat));
    AddHitable(new XZRect(0, 555, 0, 555, 0, white_mat));
    AddHitable(new FlipFace(new XZRect(0, 555, 0, 555, 555, white_mat)));
    AddHitable(new FlipFace(new XYRect(0, 555, 0, 555, 555, white_mat)));

    AddHitable(new Box(Vec3(130, 0, 65), Vec3(295, 165, 230), white_mat));
    AddHitable(new Box(Vec3(265, 0, 295), Vec3(430, 330, 460), white_mat));
  }
};

#endif  // __SCENE_CORNELL_BOX_H__
