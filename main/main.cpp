#include <time.h>
#include <stdlib.h>
#include <iostream>

#include "RandomWorld.h"

static void run_scene() {
    auto scene = new TestRandomWorld();
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