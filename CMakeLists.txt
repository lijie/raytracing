# cmake

project(raytracing)

# add_executable(rtapp main/ray_tracing.cpp)
# target_include_directories(rtapp PRIVATE src)

add_executable(rtmain main/main.cpp src/scene.cpp src/render.cpp src/rand.cpp src/pdf.cpp)
target_include_directories(rtmain PRIVATE src scenes)

add_executable(mctest main/mctest.cpp)
