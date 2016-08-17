#include <iostream>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

#include "msg.h"
#include "cosmology.h"
#include "particle.h"
#include "halo.h"
#include "satellite.h"
#include "rand.h"

using namespace std;


static const int ninterp= 1001;
static const double xmax= 100;
static double spline_ymax= 0.0;
static gsl_spline* spline= 0;
static gsl_interp_accel* acc= 0;

static void init_I();
static gsl_spline* spline_I= 0;
static gsl_interp_accel* acc_I= 0;
static double xmin_I, ymin_I, xmax_I, ymax_I;


static double f_inverse(const double fx);
static double compute_v_rms(const double r,
		       const double Mvir, const double rvir, const double cvir);
static void random_direction(float e[]);

static inline double f(const double x)
{
  return log(1.0 + x) - x/(1.0 + x);
}

void satellite_init()
{
  msg_printf(msg_debug, "start satellite_init()\n");
  
  double* x= (double*) malloc(sizeof(double)*2*ninterp); assert(x);
  double* y= x + ninterp;
  for(int i=0; i<ninterp; ++i) {
    x[i]= xmax*i/(ninterp-1);
    y[i]= f(x[i]);
  }

  spline_ymax= y[ninterp-1];

  spline= gsl_spline_alloc(gsl_interp_cspline, ninterp);
  acc= gsl_interp_accel_alloc(); 
  gsl_spline_init(spline, y, x, ninterp);

  init_I();
  
  free(x);
  msg_printf(msg_verbose, "satellite module initialised\n");
}

void satellite_free()
{
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
}

void satellite(Halo const * const h, Particle* const g)
{
#ifdef DEBUG
  assert(spline);
#endif

  const double a= 1.0/(1.0 + h->z);
  const double rho_m= cosmology_rho_m()/(a*a*a); // physical [1/h Mpc]^-3
  const double r200m= 1000.0*pow(h->M/(4.0*M_PI/3.0*200.0*rho_m), 1.0/3.0);
    // physical 1/h kpc
  const double rs= h->rs; // physical 1/h kpc.
  const double c200m= r200m/rs;
  
  // draw random mass M(r)/M0 between [0, f(c200m)]
  const double fmax= f(c200m);
  const double fx= fmax*rand_uniform();

  // solve for f(x) = fx, where x= r/r_s
  double x= f_inverse(fx);

  double r_sat= x*rs; // location of the satellite from center [phys /h kpc]
  
  // compute vrms(r)
  double vrms= compute_v_rms(r_sat, h->M, r200m, c200m);

  r_sat= r_sat/(1000.0f*a); // physical /h kpc -> comoving /h Mpc

  float e[3];
  random_direction(e);

  // satellite x v contains only offset from halo
  
  g->x[0] = r_sat*e[0];
  g->x[1] = r_sat*e[1];
  g->x[2] = r_sat*e[2];

  g->vr= vrms*rand_gaussian();
}

//
// NFW f_inverse function
//
double f_inverse(const double fx)
{
#ifdef DEBUG
  assert(spline);

  if(!(0 <= fx && fx <= spline_ymax)) {
    fprintf(stderr, "Error: spline f_inverse out of range\n"
	    "%e spline_ymax %e\n", fx, spline_ymax);
  }
  #endif

  assert(0 <= fx && fx <= spline_ymax);
  
  double x= gsl_spline_eval(spline, fx, acc);
#ifdef DEBUG
  assert(fabs(fx - f(x)) < 0.1);
#endif
  
  return x;
}


//
// NFW velocity rms
//   See, e.g. de la Torre et al. 2013 (A.6)

static double g(const double x, void* param)
{
  return f(x)/(x*x*x*(1.0+x)*(1.0+x));
}
    

double compute_v_rms(const double r,
		     const double Mvir, const double rvir, const double cvir)
{
  // ToDo: interpolation would be faster for I.

  
  const double G=43007.1/1.0e10; // for 1/h kpc, solar mass, km/s
  const double rs= rvir/cvir;
  const double s1= 1.0 + r/rs;


  double x = r/rs;
  double I= ymin_I;
  if(x > xmin_I)
    I= gsl_spline_eval(spline_I, r/rs, acc_I);
  return sqrt(G*Mvir/(f(cvir)*rs)*(r/rs)*s1*s1*I);
  //return sqrt(G*Mvir/(f(cvir)*rs)*(r/rs)*s1*s1*I(r/rs));
}

void init_I()
{
  const int worksize= 1000;
  gsl_integration_workspace *w= gsl_integration_workspace_alloc(worksize);
      
  gsl_function F;
  F.function = &g;
  F.params = 0;

  xmin_I= 0.01;
  xmax_I= 100.0;
  const int ninterp= 200;

  double* const x= (double*) malloc(sizeof(double)*ninterp*2);
  double* const y= x + ninterp;

  spline_I= gsl_spline_alloc(gsl_interp_cspline, ninterp);
  acc_I= gsl_interp_accel_alloc(); 

  const double log_x_min= log(xmin_I);
  const double log_x_max= log(xmax_I);

  double result, abserr;
  
  for(int i=0; i<ninterp; ++i) {
    x[i]= exp(log_x_min + (log_x_max - log_x_min)*i/(ninterp-1));


    gsl_integration_qagiu(&F, x[i], 0, 1.0e-5, worksize, w, &result, &abserr);
    y[i]= result;
  }

  ymin_I= y[0];
  
  gsl_spline_init(spline_I, x, y, ninterp);

  gsl_integration_workspace_free(w);
  free(x);
}


//
// random direction
//
void random_direction(float e[])
{
  float y1, y2, r2;

  do{
    y1= 1.0f-2.0f*rand_uniform();
    y2= 1.0f-2.0f*rand_uniform();
    r2= y1*y1 + y2*y2;
  } while(r2 > 1.0f);

  float sq1r= sqrt(1.0f-r2);
  e[0]= 2*y1*sq1r;
  e[1]= 2*y2*sq1r;
  e[2]= 1.0f - 2.0f*r2;
}
