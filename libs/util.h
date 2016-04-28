#ifndef UTIL_H
#define UTIL_H 1

#include <cmath>

namespace util {
  inline float norm(const float x[]) {
    return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
  }
}

#endif
