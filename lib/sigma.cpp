// sigma(M) is rms fluctuation of density fluctuation smoothed on scale mass M

#include <iostream>
#include <cmath>
#include <cassert>

#include <gsl/gsl_interp.h>

#include "const.h"
#include "sigma.h"
#include "power.h"
#include "cosmology.h" // => rho_m

using namespace std;

static int n= 0;
static double *M= 0;
static double *sinv;
static double M_min, M_max;
static double sinv_min, sinv_max;
static gsl_interp *interp, *interp_inv;
static gsl_interp_accel *acc, *acc_inv;


void sigma_init(const double M_min_, const double M_max_, const int n_)
{
  if(M) return;
  
  const double rho_m= cosmology_rho_m(); assert(rho_m >= 0.0);

  M_min= M_min_;
  M_max= M_max_;
  n= n_;
  const double log_Mmin= 0.99*log(M_min);
  const double log_Mmax= 1.01*log(M_max);

  M= (double*) malloc(sizeof(double)*2*n); assert(M);
  sinv= M + n;
  
  for(int i=0; i<n; ++i) {
    double logM = log_Mmin + (log_Mmax - log_Mmin)*((double) i/(n-1));
    M[i]= exp(logM);
    double R= pow(M[i]/(4.0*M_PI/3.0*rho_m), 1.0/3.0);

    sinv[i]= 1.0/power_compute_sigma(R);
  }

  // Function: 1/sigma -> M
  interp = gsl_interp_alloc(gsl_interp_cspline, n);
  gsl_interp_init(interp, sinv, M, n);
  acc = gsl_interp_accel_alloc(); 

  // Function: M -> 1/sigma
  interp_inv = gsl_interp_alloc(gsl_interp_cspline, n);
  gsl_interp_init(interp_inv, M, sinv, n); 
  acc_inv = gsl_interp_accel_alloc();

  sinv_min= sigma_inv(M_min);
  sinv_max= sigma_inv(M_max);
}

void sigma_free()
{
  gsl_interp_free(interp);
  gsl_interp_free(interp_inv);
  gsl_interp_accel_free(acc);
  gsl_interp_accel_free(acc_inv);

  free(M);
}

bool sigma_initilised()
{
  return !(M == 0);
}

double sigma_M(const double sigma0)
{
  // sigma0 is sigma(M, z=0)
#ifdef DEBUG
  assert(n > 0);
  if(1/sigma0 < sinv_min || 1/sigma0 > sinv_max) {
    cerr << "Error: sigma_M(sigma0) out of range.\n";
    cerr << sinv[0] << " " << sinv[n-1] << " " << 1.0/sigma0 << endl;
  }
#endif

  return gsl_interp_eval(interp, sinv, M, 1.0/sigma0, acc);
}

double sigma_inv(const double MM)
{
#ifdef DEBUG
  assert(n > 0);
  assert(M_min <= MM && MM <= M_max);
#endif
  
  return gsl_interp_eval(interp_inv, M, sinv, MM, acc_inv);
}

double* sigma_M_array()
{
  assert(M);
  return M;
}

double* sigma_sinv_array()
{
  assert(M);
  return sinv;
}

double sigma_M_min()
{
  assert(M);
  return M_min;
}

double sigma_M_max()
{
  assert(M);
  return M_max;
}

double sigma_sinv_min()
{
  assert(M);
  return sinv_min;
}

double sigma_sinv_max()
{
  assert(M);
  return sinv_max;
}

int sigma_n()
{
  return n;
}
