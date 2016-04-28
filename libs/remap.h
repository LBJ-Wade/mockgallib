#ifndef REMAP_H
#define REMAP_H 1

#include <vector>
#include "halo.h"


struct Remapping {
  float e[9]; // orthonormal basis of remapped coordinate system
  float l[3]; // boxsize[3]/boxsize_original
  float boxsize[3];
  int icopy_begin[3], icopy_end[3];
};

//void remap_basis(const int u[], float* const e, float* const l);

Remapping* remap_alloc(const int u[], const float boxsize);
void remap_free(Remapping* const remap);

//void remap_coordinate(std::vector<Halo>& v, const float e[], const float l[], const float boxsize);
void remap_get_boxsize(Remapping const * const remap, float * const boxsize);
void remap_halo_coordinate(Halo* const h, Remapping const * const remap, const float boxsize);
//void remap_coordinate(std::vector<Halo>& v, Remapping const * const remap, const float boxsize);

class ErrorRemap {
};


#endif
