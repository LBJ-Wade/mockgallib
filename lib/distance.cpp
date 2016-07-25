//
// Compute redshift from distance r
//
// Dependence: cosmology::omega_m
//
#include <cstdlib>
#include <cassert>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>

#include "msg.h"
#include "cosmology.h"
#include "distance.h"


using namespace std;

static gsl_interp_accel *acc= 0;
static gsl_spline* spline= 0;
static double z_max, d_max= 0.0;

void distance_init(const double z_max_)
{
  if(acc) {
    if(z_max == z_max_)
      return;
    else
      msg_abort("distance already initialzed");
  }

  static size_t n= 1001;

  spline= gsl_spline_alloc(gsl_interp_cspline, n);
  acc= gsl_interp_accel_alloc();  
  z_max= z_max_;

  double* const d= (double*) malloc(sizeof(double)*n*2); assert(d);
  double* const z= d + n;

  d[0]= 0.0; z[0]= 0.0;

  for(int i=1; i<n; ++i) {
    z[i]= z_max*i/(n-1);
    double a= 1.0/(1.0 + z[i]);
    
    d[i]= cosmology_compute_comoving_distance(a);
  }

  d_max= d[n-1];
  
  gsl_spline_init(spline, d, z, n);
  msg_printf(msg_debug, "distance_init max distance %.2f\n", d[n-1]);

  free(d);
}

void distance_free()
{
  if(spline) gsl_spline_free(spline);
  if(acc) gsl_interp_accel_free(acc);
}

double distance_redshift(const double d)
{
#ifdef DEBUG
  assert(spline);
#endif
  if(d > d_max) return -1.0;
  
  return gsl_spline_eval(spline, d, acc);
}


double distance_max()
{
  return d_max;
}
