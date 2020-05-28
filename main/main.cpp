#include <time.h>
#include <stdlib.h>
#include <iostream>

#include "RandomWorld.h"
#include "RectLight.h"

static void run_scene() {
    // auto scene = new TestRandomWorld();
    auto scene = new RectLightScene();
    scene->Create();
    scene->Run();
}

int main() {
  srand(time(NULL));
  std::cout << time(NULL) << std::endl;
  run_scene();
  std::cout << time(NULL) << std::endl;
  return 0;
}