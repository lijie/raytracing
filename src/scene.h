#ifndef __SCENE_H__
#define __SCENE_H__

#include <vector>

#include "ray.h"

class Hitable;
class Camera;

class Scene {
 public:
  bool Create();
  void Destroy();
  void Run();

  void AddHitable(Hitable *obj);
  Hitable *GetRoot() const { return world_; }
  void SetOutput(const char *path);
  void SetCamera(Camera *c) { camera_ = c; };
  void SetScreenSize(double width, double height) {
    width_ = width;
    height_ = height;
  }
  void SetLightSource(Hitable *light) { light_ = light; }

  virtual Vec3 BackgroundColor(const Ray &ray) const;
  Hitable *GetLightSource() const { return light_; }

 protected:
  virtual void OnCreate() {}
  virtual void OnDestroy() {}

 private:
  Hitable *world_;
  Hitable *light_;
  Camera *camera_;
  double width_;
  double height_;
  std::vector<Hitable *> pending_vec_;
  const char *output_file_path_;
};

#endif  // __SCENE_H__
