#ifndef SKY_H
#define SKY_H 1

#include <vector>
#include <cassert>
#include "halo.h"
#include "util.h"

class Sky {
 public:
  Sky(const double ra[], const double dec[], const double z[] );
  inline void compute_radec(const float x[], float* const radec) const;
  double ra_range[2], dec_range[2], r_range[2];
  double r_min, r_max;
  float left[3], right[3], width[3]; // bounding box
  float ra0, dec0, theta0, cos_theta0, sin_theta0;
};

void Sky::compute_radec(const float x[], float* const radec) const
{
  // compute h->ra and dec from cuboid coordinate h->x
  // ra0, dec0: RA and Dec of xaxis (y=0, z=0)

  float theta= asin(x[2]/util::norm(x));
    
  // rotate cosO in x-z plane
  float x1= cos_theta0*x[0] - sin_theta0*x[2];
  float y1= x[1];
  float rr= sqrt(x1*x1 + y1*y1);
  float phi= asin(y1/rr); if(x1 < 0) phi= M_PI - phi;

  radec[0]= ra0  - 180.0/M_PI*phi; // RA is right to left
  radec[1]= dec0 + 180.0/M_PI*theta;
    
  assert(-180.0f <= radec[1] && radec[1] <= 180.0f);
}


  
#endif
