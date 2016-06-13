#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>

#include "msg.h"
#include "remap.h"
#include "util.h"
// Remapping of a cubic box to a cuboid
// Carlson & White, 2010, ApJS, 190, 311

using namespace std;

static int no_remap= 0;

static inline void get_icopy(const float e[], const float l[],
	       int* const icopy_begin, int* const icopy_end);

static inline void rotate_vector(float * const vec, const float e[])
{
  float v[3];
  v[0]= util::dot(vec, e);
  v[1]= util::dot(vec, e+3);
  v[2]= util::dot(vec, e+6);

  vec[0]= v[0];
  vec[1]= v[1];
  vec[2]= v[2];
}

static inline void assert_float(const float x, const float x0,
				const float eps=1.0e-4)
{
  assert(fabs(x - x0) < eps);
}

//
// Returns orthornormal basis e[] and boxsize after remapping l[]
// for given integer vectors u[9];
// 3x3 matrix u must have det=1
//

Remap::Remap(const int u[], const float cboxsize)
{
  cboxsize_= cboxsize;
  
  int det = u[0]*u[4]*u[8] - u[0]*u[5]*u[7] - u[1]*u[3]*u[8]
              + u[1]*u[5]*u[6] + u[2]*u[3]*u[7] - u[2]*u[6]*u[4];
  if(det != 1) {
    msg_printf(msg_fatal,
	       "Error: determinant of u matrix must be 1. det= %d\n"
	       "u= %d %d %d %d %d %d %d %d %d\n",
	       det, u[0], u[1], u[2], u[3], u[4], u[5], u[6], u[7], u[8]);
    throw ErrorRemap();
  }

  // orthognal vectors (before normalisation)
  float w[9];

  // w_1 = u_1
  w[0]= (float) u[0];
  w[1]= (float) u[1];
  w[2]= (float) u[2];

  // w_2 = u_2 - u_1.u_2/|u_1|^2

  const float alpha= -float(u[0]*u[3] + u[1]*u[4] + u[2]*u[5])/
                      float(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);

  w[3]= u[3] + alpha*u[0];
  w[4]= u[4] + alpha*u[1];
  w[5]= u[5] + alpha*u[2];

  // w_3 = w_1 * w_2
  w[6]= w[1]*w[5] - w[2]*w[4];
  w[7]= w[2]*w[3] - w[0]*w[5];
  w[8]= w[0]*w[4] - w[1]*w[3];

  const float n1= util::norm(w);
  const float n2= util::norm(w+3);
  const float n3= util::norm(w+6);

  // cuboid boxsize in units of the original boxsize
  l_[0]= n1;
  l_[1]= n2;
  l_[2]= 1.0/(n1*n2);

  for(int k=0; k<3; ++k)
    this->boxsize[k]= cboxsize*l_[k];

  // orthonormal basis
  float* const e= this->e_;

  e[0]= w[0]/n1;
  e[1]= w[1]/n1;
  e[2]= w[2]/n1;
  e[3]= w[3]/n2;
  e[4]= w[4]/n2;
  e[5]= w[5]/n2;
  e[6]= w[6]/n3;
  e[7]= w[7]/n3;
  e[8]= w[8]/n3;
  
  assert_float(util::dot(e, e), 1.0f);
  assert_float(util::dot(e, e+3), 0.0f);
  assert_float(util::dot(e, e+6), 0.0f);
  assert_float(util::dot(e+3, e+3), 1.0f);
  assert_float(util::dot(e+3, e+6), 0.0f);
  assert_float(util::dot(e+6, e+6), 1.0f);

  // Find necessary range of periodic copies of the simulation cube
  get_icopy(e, l_, icopy_begin_, icopy_end_);
}

bool Remap::coordinate(Halo* const h) const
{
  //
  // Rotate the halo coordinate to a cuboid
  // Input:
  //   h:    Halo
  //   e[9]: orthonormal basis of cuboid
  //   l[3]: cuboid boxsize / boxsize
  //   boxsize: boxsize of the original cube
  //
  // Output:
  //   h->x[3] is converted to cuboid coordinate
  //

  float const * const e= e_;
  float const * const l= l_;

  // find necessary periodic copies of the cube

  const float lx= l[0]*cboxsize_;
  const float ly= l[1]*cboxsize_;
  const float lz= l[2]*cboxsize_;

  float x[3], r[3];

  for(int ix=icopy_begin_[0]; ix<icopy_end_[0]; ++ix) {
    x[0]= h->x[0] + ix*cboxsize_;
    for(int iy=icopy_begin_[1]; iy<icopy_end_[1]; ++iy) {
      x[1]= h->x[1] + iy*cboxsize_;
      for(int iz=icopy_begin_[2]; iz<icopy_end_[2]; ++iz) {
	x[2]= h->x[2] + iz*cboxsize_;

	r[0]= util::dot(x, e);
	r[1]= util::dot(x, e+3);
	r[2]= util::dot(x, e+6);

	if(0.0f <= r[0] && r[0] < lx &&
	   0.0f <= r[1] && r[1] < ly &&
	   0.0f <= r[2] && r[2] < lz) {

	  h->x[0]= r[0]; h->x[1]= r[1]; h->x[2]= r[2];
	  
	  rotate_vector(h->v, e);
	  //rotate_vector(p->f, e);

	  return true;
	}
      }
    }
  }
  no_remap++;

  // remapping not found
  h->x[0]= 0; h->x[1]= 0; h->x[2]= 0;

  return false;
}


void get_icopy(const float e[], const float l[],
	       int* const icopy_begin, int* const icopy_end)
{
  // return range of periodic replication necessary to cover the remapped
  // volume
  
  for(int k=0; k<3; ++k) {
    // corners of cuboid
    float xc[8];
    xc[0]= 0.0f;
    xc[1]= xc[0] + l[0]*e[0+k];
    xc[2]= xc[0] + l[1]*e[3+k];
    xc[3]= xc[0] + l[2]*e[6+k];
    xc[4]= xc[0] + l[0]*e[0+k] + l[1]*e[3+k];
    xc[5]= xc[0] + l[0]*e[0+k] + l[2]*e[6+k];
    xc[6]= xc[0] + l[1]*e[3+k] + l[2]*e[6+k];
    xc[7]= xc[0] + l[0]*e[0+k] + l[1]*e[3+k] + l[2]*e[6+k];
    
    float xmin= xc[0], xmax= xc[0];
    for(int i=1; i<7; ++i) {
      if(xc[i] < xmin) xmin= xc[i];
      if(xc[i] > xmax) xmax= xc[i];
    }

    icopy_begin[k] = (int)floor(xmin);
    icopy_end[k] = (int)ceil(xmax);
  }
}

