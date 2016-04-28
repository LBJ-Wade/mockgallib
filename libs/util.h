#ifndef UTIL_H
#define UTIL_H 1

#include <cmath>

namespace util {
  static inline float norm(const float x[]) {
    return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
  }

  static inline float dot(float const * const x, float const * const y) {
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
  }

}

#endif
