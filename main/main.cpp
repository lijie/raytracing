#include <time.h>
#include <stdlib.h>
#include <iostream>

#include "RandomWorld.h"
#include "RectLight.h"
#include "CornellBox.h"

static void run_scene() {
    // auto scene = new TestRandomWorld();
    // auto scene = new RectLightScene();
    auto scene = new CornellBox();
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