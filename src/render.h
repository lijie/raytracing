#ifndef __RENDER_H__
#define __RENDER_H__

#include <vector>

#include "vec3.h"

class Scene;
class Camera;
class Ray;

class Renderer {
 private:
  int nr_of_thread_ = 8;

 public:
  Renderer(int nr_of_thread = 8): nr_of_thread_(nr_of_thread) {}

  // depth 控制反射次数
  virtual Vec3 color_of_ray(const Ray& ray, const Scene *scene, int depth, bool trace = false);

  // Render entire world
  void Render(int nx, int ny, const Camera& camera, const Scene *scene,
              std::vector<int>* data);

};
#endif
