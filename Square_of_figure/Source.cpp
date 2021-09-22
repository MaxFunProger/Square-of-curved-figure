#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <string.h>
#include <stack>
#include <queue>
#include <cmath>
#include <random>

const double kEps = 1e-7;
const double kEps_c = 1e-5;

double first_func(double x) {
  return 3 / ((x - 1) * (x - 1) + 1);
}

double first_func_d(double x) {
  return -(6 * (-1 + x)) / ((1 + (-1 + x) * (-1 + x)) * (1 + (-1 + x) * (-1 + x)));
}

double first_func_d2(double x) {
  return (6 * (2 - 6 * x + 3 * x * x)) / ((2 - 2 * x + x * x) * (2 - 2 * x + x * x) * (2 - 2 * x + x * x));
}

double second_func(double x) {
  return sqrt(x + 0.5);
}

double second_func_d(double x) {
  return 1 / (2 * sqrt(x + 0.5));
}

double second_func_d2(double x) {
  return -(x + 0.5) * (x + 0.5) / 4 / (x + 0.5) / (x + 0.5) / (x + 0.5);
}

double third_func(double x) {
  return exp(-x);
}

double third_func_d(double x) {
  return -exp(-x);
}

double third_func_d2(double x) {
  return exp(-x);
}

double second_third(double x) {
  return second_func(x) - third_func(x);
}

double first_second(double x) {
  return first_func(x) - second_func(x);
}

double first_third(double x) {
  return first_func(x) - third_func(x);
}

double second_third_d(double x) {
  return second_func_d(x) - third_func_d(x);
}

double first_second_d(double x) {
  return first_func_d(x) - second_func_d(x);
}

double first_third_d(double x) {
  return first_func_d(x) - third_func_d(x);
}

double second_third_d2(double x) {
  return second_func_d2(x) - third_func_d2(x);
}

double first_second_d2(double x) {
  return first_func_d2(x) - second_func_d2(x);
}

double first_third_d2(double x) {
  return first_func_d2(x) - third_func_d2(x);
}

double roots(double a, double b, double (*func)(double), double (*func_d)(double), double(*func_d2)(double)) {
  while (std::fabs(a - b) > kEps) {
    if (func(a) * func_d2(a) < 0) {
      a = a - func(a) * (a - b) / (func(a) - func(b));
    }
    else {
      a = a - func(a) / func_d(a);
    }

    if (func(b) * func_d2(b) < 0) {
      b = b - func(b) * (b - a) / (func(b) - func(a));
    }
    else {
      b = b - func(b) / func_d(b);
    }
      
  }
  return (a + b) / 2;
}


double square(double a, double b, double (*func)(double)) {
  int n = 1024;
  double step = (b - a) / n;
  double s1 = 0;
  double s2 = 2 * kEps_c;
  while (std::fabs(s2 - s1) > kEps_c) {
    step = (b - a) / n;
    s1 = s2;
    s2 = 0;
    double cnt = a + step;
    double pre = a;
    while (cnt <= b) {
      s2 += (cnt - pre) / 6 * (func(pre) + 4 * func((pre + cnt) / 2) + func(cnt));
      pre = cnt;
      cnt += step;
    }
    if (pre < b && cnt > b) {
      s2 += (b - pre) / 6 * (func(pre) + 4 * func((pre + b) / 2) + func(b));
    }
    n *= 2;
  }
  return s2;
}

double fRand(double fMin, double fMax)
{
  double f = (double)rand() / RAND_MAX;
  return fMin + f * (fMax - fMin);
}



double casino(double x1, double x2, double x3, double (*func1)(double), double (*func2)(double), double (*func3)(double)) {
  int all = 0;
  int inside = 0;
  for (int i = 0; i < 1'000'000; ++i) {
    double x = fRand(x1, x3);
    double y = fRand(0, 3);
    ++all;
    if (y <= func1(x) && y >= func2(x) && y >= func3(x))
      ++inside;
  }
  return (x3 - x1) * 3 * inside / all;


}






int main() {
  double x1 = roots(-1, 0, first_third, first_third_d, first_third_d2);
  double x2 = roots(0, 1, second_third, second_third_d, second_third_d2);
  double x3 = roots(1, 3, first_second, first_second_d, first_second_d2);
  std::cout.precision(10);
  std::cout << "Roots:" << "\n";
  std::cout << "\t" << x1 << "\n";
  std::cout << "\t" << x2 << "\n";
  std::cout << "\t" << x3 << "\n";
  double result = square(x1, x3, first_func) - square(x1, x2, third_func) - square(x2, x3, second_func);
  std::cout << "Simpson:" << "\n";
  std::cout << "\t" << result << "\n";
  std::cout << "Casino:" << "\n";
  std::cout << "\t" << casino(x1, x2, x3, first_func, second_func, third_func);


  return 0;
}