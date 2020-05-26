#include "scene.h"

#include <fstream>
#include <iostream>

#include "render.h"

static void write_ppm(const char* filename, int nx, int ny,
                      const std::vector<int>& data) {
  std::ofstream out;

  out.open(filename);

  out << "P3\n" << nx << " " << ny << "\n255\n";  // PPM header

  for (size_t i = 0; i < data.size(); i++) {
    out << data[i] << "\n";
  }

  out.close();
}

bool Scene::Create() {
  OnCreate();
  return true;
}

void Scene::Destroy() { OnDestroy(); }

void Scene::SetOutput(const char* path) { output_file_path_ = path; }

void Scene::AddHitable(Hitable* obj) { pending_vec_.push_back(obj); }

void Scene::Run() {
  Renderer renderer;
  Hitable* world = GetRoot();
  std::vector<int> data;
  renderer.Render(width_, height_, *camera_, world, &data);
  write_ppm(output_file_path_, width_, height_, data);
  Destroy();
}