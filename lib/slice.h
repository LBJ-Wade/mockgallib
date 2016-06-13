#ifndef SLICE_H
#define SLICE_H 1

#include <vector>
#include <cmath>
#include "halo.h"

class Slice {
 public:
  Slice(const float boxsize[], const float sky_box[], const float center[]);
  void transform(Halo* const h) const;

  float boxsize[3];
  int n;
 private:
  int n_[3];     // number of subvolumes that fits the cuboid for each direction
  int coord_[3]; // x[coord[i]]
  float x0_[3];  // left cornor coordinate of position after slicing
};

#endif
