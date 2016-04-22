#include <cstdio>
#include <cmath>
#include <cassert>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>

#include "msg.h"
#include "cosmology.h"
#include "particle.h"
#include "halo.h"
#include "satellite.h"

static gsl_rng* random_generator= 0;
static gsl_root_fdfsolver* solver= 0;

static double f_inverse(const double fx, double x);
static double compute_v_rms(const double r,
		       const double Mvir, const double rvir, const double cvir);
static void random_direction(float e[]);

static inline double f(const double x)
{
  return log(1.0+x) - x/(1.0+x);
}

void satellite_init(gsl_rng* const rng)
{
  random_generator= rng;

  if(!solver)
    solver= gsl_root_fdfsolver_alloc(gsl_root_fdfsolver_newton);
}

void satellite_free()
{
  if(!solver) {
    gsl_root_fdfsolver_free(solver);
    solver= 0;
  }
}

void satellite(Halo const * const h, Particle* const g)
{
#ifdef DEBUG
  assert(random_generator);
  assert(solver);
#endif
  
  const double a= 1.0/(1.0 + h->z);
  const double rho_m= cosmology_rho_m()/(a*a*a); // physical [1/h Mpc]^-3
  const double r200m= 1000.0*pow(h->M/(4.0*M_PI/3.0*200.0*rho_m), 1.0/3.0);
    // physical 1/h kpc
  const double c200m= r200m/h->rs;

  //fprintf(stderr, "r200m c rs %e %e %e\n", r200m, c200m, h->rs);

  // draw random mass M(r)/M0 between [0, f(c200m)]
  const double fmax= f(c200m);
  const double fx= fmax*gsl_rng_uniform(random_generator);

  // solve for f(x) = fx, where x= r/r_s
  double x= c200m*fx/fmax; // initial guess

  x= f_inverse(fx, x);

  double r_sat= x*h->rs; // location of the satellite from center
  
  // compute vrms(r)
  double vrms= compute_v_rms(r_sat, h->M, r200m, c200m);

  r_sat= r_sat/(1000.0f*a); // physical /h kpc -> comoving /h Mpc

  float e[3];
  random_direction(e);

  // satellite x v contains only offset from halo
  
  g->x[0] = r_sat*e[0];
  g->x[1] = r_sat*e[1];
  g->x[2] = r_sat*e[2];

  g->vr= vrms*gsl_ran_ugaussian(random_generator);
}

//
// NFW f_inverse function
//
double solver_f(double x, void* fx)
{
  return f(x) - *(double*)fx;
}

double solver_df(double x, void* fx)
{
  double x1= x + 1.0;
  return x/(x1*x1);
}

void solver_fdf(double x, void* fx, double* f_out, double* df_out)
{
  double x1= x + 1.0;
  *f_out= f(x) - *(double*)fx;
  *df_out= x/(x1*x1);
}

double f_inverse(const double fx, double x)
{
  // x is the initial guess of the root f(x)=fx
  gsl_function_fdf fdf;
  fdf.f= &solver_f;
  fdf.df= &solver_df;
  fdf.fdf= &solver_fdf;
  fdf.params= (void*) &fx;

  gsl_root_fdfsolver_set(solver, &fdf, x);

  int iter= 0, status;
  double x0;
  do {
    iter++;
    status= gsl_root_fdfsolver_iterate(solver);
    x0= x;
    x= gsl_root_fdfsolver_root(solver);
    status= gsl_root_test_delta(x, x0, 0, 1.0e-6);

    if(status == GSL_SUCCESS)
      break;
  } while(status == GSL_CONTINUE && iter < 100);

  double fx_check= f(x);

  assert(fabs(fx - fx_check) < 1.0e-4);
  
  return x;
}


//
// NFW velocity rms
//   See, e.g. de la Torre et al. 2013 (A.6)

static double g(const double x, void* param)
{
  return f(x)/(x*x*x*(1.0+x)*(1.0+x));
}
    

static double I(const double x)
{
  const int worksize= 1000;
  gsl_integration_workspace *w= gsl_integration_workspace_alloc(worksize);

  gsl_function F;
  F.function = &g;
  F.params = 0;

  double result, abserr;
  gsl_integration_qagiu(&F, x, 0, 1.0e-8, worksize, w, &result, &abserr);

  return result;
}

double compute_v_rms(const double r,
		     const double Mvir, const double rvir, const double cvir)
{
  const double G=43007.1/1.0e10; // for 1/h kpc, solar mass, km/s
  const double rs= rvir/cvir;
  const double s1= 1.0 + r/rs;
  return sqrt(G*Mvir/(f(cvir)*rs)*(r/rs)*s1*s1*I(r/rs));
}

//
// random direction
//
void random_direction(float e[])
{
  float y1, y2, r2;

  do{
    y1= 1.0f-2.0f*gsl_rng_uniform(random_generator);
    y2= 1.0f-2.0f*gsl_rng_uniform(random_generator);
    r2= y1*y1 + y2*y2;
  } while(r2 > 1.0f);

  float sq1r= sqrt(1.0f-r2);
  e[0]= 2*y1*sq1r;
  e[1]= 2*y2*sq1r;
  e[2]= 1.0f-2.0f*r2;
}
