#include <stdio.h>

#include <random>

// 推算 Pi
int estimate_pi_1() {
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  int N = 100000000;
  int in_circle = 0;

  for (int i = 0; i < N; i++) {
    auto x = dist(mt);
    auto y = dist(mt);

    if (x * x + y * y < 1) {
      in_circle++;
    }
  }

  printf("Pi is %f\n", 4 * ((double)in_circle / (double)N));
  return 0;
}

// 单纯的随机随着N增大, 存在收益递减问题
// 把整个估计面积划分为 sqrt_N * sqrt_N 个格子
// 在小格子内随机, 使得我们的随机数能均匀分布到整个估计面积
// 同样的循环次数, 比单纯的随机有着更高的收益
int estimate_pi_2() {
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dist;
  int sqrt_N = 10000;
  int in_circle = 0;

  for (int i = 0; i < sqrt_N; i++) {
    for (int j = 0; j < sqrt_N; j++) {
      auto x = (i + dist(mt)) / sqrt_N;
      auto y = (j + dist(mt)) / sqrt_N;

      if (x * x + y * y < 1) {
        in_circle++;
      }
    }
  }

  printf("Pi is %f\n", 4 * ((double)in_circle / ((double)sqrt_N * sqrt_N)));
  return 0;
}

int main() {
  estimate_pi_1();
  estimate_pi_2();
}
