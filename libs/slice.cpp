#include <cstdlib>
#include <cmath>
#include <cassert>
#include "slice.h"

using namespace std;

Slice::Slice(const float boxsize[], const float sky_box[], const float center[] )
{
  int n1[3], n2[3];

  // Find the number of sky_boxes that fits boxsize
  
  // x y z
  n1[1]= (int) floor(boxsize[1]/sky_box[1]);
  n1[2]= (int) floor(boxsize[2]/sky_box[2]);

  // x z y
  n2[1]= (int) floor(boxsize[2]/sky_box[1]);
  n2[2]= (int) floor(boxsize[1]/sky_box[2]);


  n_[0]= (int) floor(boxsize_[0]/sky_box[0]);
  
  if(n2[1]*n2[2] > n1[1]*n1[2]) {
    // We can fit more slices by flipping y and z coordinates
    coord_[0]= 0; coord_[1]= 2; coord_[2]= 1;
    n_[1]= n2[1]; n_[2]= n2[2];
  }
  else {
    coord_[0]= 0; coord_[1]= 1; coord_[2]= 2;
    n_[1]= n1[1]; n_[2]= n1[2];
  }

  for(int k=0; k<3; ++k) {
    boxsize_[k]= boxsize[k]/n_[coord_[k]];
    x0_[k]= center[k] - 0.5f*boxsize_[k];
  }
}

