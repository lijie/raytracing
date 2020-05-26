#ifndef __SCENE_H__
#define __SCENE_H__

#include <vector>

class Hitable;
class Camera;

class Scene {
 public:
  bool Create();
  void Destroy();
  void Run();

  void AddHitable(Hitable *obj);
  Hitable *GetRoot();
  void SetOutput(const char *path);
  void SetCamera(Camera *c) { camera_ = c; };
  void SetScreenSize(double width, double height) {
    width_ = width;
    height_ = height;
  }

 protected:
  virtual void OnCreate() {}
  virtual void OnDestroy() {}

 private:
  Hitable *world_;
  Camera *camera_;
  double width_;
  double height_;
  std::vector<Hitable *> pending_vec_;
  const char *output_file_path_;
};

#endif  // __SCENE_H__
