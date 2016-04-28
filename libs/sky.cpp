//
// Sky contains a volume information of ra/dec/r in the sky
//

#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>

#include "msg.h"
#include "cosmology.h"
#include "sky.h"


using namespace std;

static inline double rad(const double deg)
{
  return deg/180.0*M_PI;
}

static void compute_bounding_box(Sky * const sky);

Sky::Sky(const double ra[], const double dec[], const double z[])
{
  ra_range[0]= ra[0]; ra_range[1]= ra[1];
  dec_range[0]= dec[0]; dec_range[1]= dec[1];

  const double a_min= 1.0/(1.0 + z[0]);
  const double a_max= 1.0/(1.0 + z[1]);
  r_range[0]= cosmology_compute_comoving_distance(a_min);
  r_range[1]= cosmology_compute_comoving_distance(a_max);

  compute_bounding_box(this);

  // The x-axis (y=0, z=0) corresponds passes the centre of RA and Dec
  dec0= dec_range[0] + 0.5*(dec_range[1] - dec_range[0]);
  ra0=  ra_range[0] + 0.5*(ra_range[1] - ra_range[0]);

  theta0= rad(dec0);
  cos_theta0= cos(theta0);
  sin_theta0= sin(theta0);
}
  



static inline void set_cartisian(const double r, const double theta, const double phi, double * const x)
{
  // Returns cartisian coordinate corresponds to polar coordinate
  // theta is measured from xy-plane, not from z axis
  const double rcosO= r*cos(theta);
  x[0]= rcosO*cos(phi);
  x[1]= rcosO*sin(phi);
  x[2]= r*sin(theta);
}

static inline void rotate_xz(double * const x, const double theta)
{
  // Rotate the vector by an angle theta (radian) in x-z plane
  const double cosO= cos(theta);
  const double sinO= sin(theta);
  
  double x1= cosO*x[0] - sinO*x[2];
  double z1= sinO*x[0] + sinO*x[2];

  x[0]= x1;
  x[2]= z1;
}

void compute_bounding_box(Sky* const sky)
{
  //
  // computes the size of cuboid that can enclose the sky region
  //
  // Output: sky->box[3]
  //
  const double theta= rad(0.5*fabs(sky->dec_range[1] - sky->dec_range[0]));
  const double phi= rad(0.5*(sky->ra_range[1] - sky->ra_range[0]));
  assert(phi> 0.0);
  assert(phi < 0.5*M_PI); // phi > pi/2 not implemented in this code

  const double theta_min= rad(sky->dec_range[0]);
  const double theta_max= rad(sky->dec_range[1]);
  const double theta0= 0.5*(theta_min + theta_max);
  double x[3];

  // x-range
  double x_min= sky->r_min;
  
  set_cartisian(sky->r_min, theta_min, phi, x);
  rotate_xz(x, -theta0);
  x_min= fmin(x_min, x[0]);
  
  set_cartisian(sky->r_min, theta_max, phi, x);
  rotate_xz(x, -theta0);
  x_min= fmin(x_min, x[0]);

  sky->right[0]= sky->r_max;
  sky->left[0]= x_min;
  sky->width[0]= sky->r_max - x_min;

  // y-range
  double y_max= 0.0;
  
  set_cartisian(sky->r_max, theta_min, phi, x);
  rotate_xz(x, -theta0);
  y_max= fmax(y_max, x[1]);
  
  set_cartisian(sky->r_max, theta_max, phi, x);
  rotate_xz(x, -theta0);
  y_max= fmax(y_max, x[1]);

  sky->left[1]= -y_max;
  sky->right[1]= y_max;
  sky->width[1]= 2.0*y_max;

  // z-range
  sky->right[2]= sky->r_max*sin(theta);
  sky->left[2]= -sky->right[2];
  sky->width[2]= 2.0*sky->right[2];
}

