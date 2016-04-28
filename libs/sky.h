#ifndef SKY_H
#define SKY_H 1

#include <vector>
#include "halo.h"
#include "util.h"

class Sky {
 public:
  Sky(const double ra[], const double dec[], const double z[] );

  double ra_range[2], dec_range[2], r_range[2];
  double r_min, r_max;
  float left[3], right[3], width[3]; // bounding box
  float ra0, dec0, theta0, cos_theta0, sin_theta0;
};




//void minimum_bounding_box(Sky const * const sky, float* const boxsize);
//void compute_ra_dec(std::vector<Halo>& v, Sky const * const sky);

template<class T>
inline void sky_compute_ra_dec(Sky const * const sky, T* const h)
{
  // compute h->ra and dec from cuboid coordinate h->x
  // ra0, dec0: RA and Dec of xaxis (y=0, z=0)

  //float const * const x= h->x;
  float theta= asin(h->x[2]/util::norm(h->x));
    
  // rotate cosO in x-z plane
  float x1= sky->cos_theta0*h->x[0] - sky->sin_theta0*h->x[2];
  float y1= h->x[1];
  float rr= sqrt(x1*x1 + y1*y1);
  float phi= asin(y1/rr); if(x1 < 0) phi= M_PI - phi;

  h->ra=  sky->ra0  - 180.0/M_PI*phi;                  // RA is right to left
  h->dec= sky->dec0 + 180.0/M_PI*theta;
    
  assert(-180.0f <= h->dec && h->dec <= 180.0f);
}


  
#endif
