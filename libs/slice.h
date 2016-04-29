#ifndef SLICE_H
#define SLICE_H 1

#include <vector>
#include <cmath>
#include "halo.h"

class Slice {
 public:
  Slice(const float boxsize[], const float sky_box[], const float center[]);
  int transform(float x[]) {
    // a copy is necessary because coord[] could change coordinates
    float xx[]= {x[0], x[1], x[2]};

    int index= (int)(((floorf(x[coord_[0]]/boxsize_[0])*n_[1]
		     + floorf(x[coord_[1]]/boxsize_[1]))*n_[0]
		     + floorf(x[coord_[2]]/boxsize_[2])));

    x[0]= x0_[0] + fmodf(xx[coord_[0]], boxsize_[0]);
    x[1]= x0_[1] + fmodf(xx[coord_[1]], boxsize_[1]);
    x[2]= x0_[2] + fmodf(xx[coord_[2]], boxsize_[2]);

    return index;
  }

  
  int n_[3];     // number of subvolumes that fits the cuboid for each direction
  int coord_[3]; // x[coord[i]]
  float boxsize_[3];
  float x0_[3];  // left cornor coordinate of position after slicing
};

#endif
