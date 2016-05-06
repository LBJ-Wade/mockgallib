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


  n_[0]= (int) floor(boxsize[0]/sky_box[0]);
  
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

  n= n_[0]*n_[1]*n_[2];
}


void Slice::transform(Halo* const h) const {
  // a copy is necessary because coord[] could change coordinates
  float x[]= {h->x[0], h->x[1], h->x[2]};
  
  int index= (int)(((floorf(x[coord_[0]]/boxsize_[0])*n_[1]
		     + floorf(x[coord_[1]]/boxsize_[1]))*n_[0]
		    + floorf(x[coord_[2]]/boxsize_[2])));
  
  h->x[0]= x0_[0] + fmodf(x[coord_[0]], boxsize_[0]);
  h->x[1]= x0_[1] + fmodf(x[coord_[1]], boxsize_[1]);
  h->x[2]= x0_[2] + fmodf(x[coord_[2]], boxsize_[2]);
  
  h->slice= index;
}

