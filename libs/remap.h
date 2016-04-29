#ifndef REMAP_H
#define REMAP_H 1

#include <vector>
#include "halo.h"


class Remap {
 public:
  Remap(const int u[], const float cboxsize);
  bool coordinate(Halo* const h) const;
  float boxsize[3];
 private:
  float cboxsize_;
  float e_[9]; // orthonormal basis of remapped coordinate system
  float l_[3]; // boxsize[3]/boxsize_original
  int   icopy_begin_[3], icopy_end_[3];
};


//void remap_get_boxsize(Remapping const * const remap, float * const boxsize);



class ErrorRemap {
};


#endif
