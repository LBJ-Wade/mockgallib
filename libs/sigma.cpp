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

Sigma::Sigma(PowerSpectrum const * const ps,
	     const double M_min1, const double M_max1, const int n1)
{
  const double rho_m= cosmology_rho_m(); assert(rho_m >= 0.0);

  M_min= M_min1;
  M_max= M_max1;
  n= n1;
  const double log_Mmin= 0.99*log(M_min);
  const double log_Mmax= 1.01*log(M_max);

  M_= (double*) malloc(sizeof(double)*2*n); assert(M_);
  sinv_= M_ + n;
  
  for(int i=0; i<n; ++i) {
    double logM = log_Mmin + (log_Mmax - log_Mmin)*((double) i/(n-1));
    M_[i]= exp(logM);
    double R= pow(M_[i]/(4.0*M_PI/3.0*rho_m), 1.0/3.0);

    sinv_[i]= 1.0/ps->compute_sigma(R);
  }

  // Function: 1/sigma -> M
  interp_ = gsl_interp_alloc(gsl_interp_cspline, n);
  gsl_interp_init(interp_, sinv_, M_, n);
  acc_ = gsl_interp_accel_alloc(); 

  // Function: M -> 1/sigma
  interp2_ = gsl_interp_alloc(gsl_interp_cspline, n);
  gsl_interp_init(interp2_, M_, sinv_, n); 
  acc2_ = gsl_interp_accel_alloc();

  sinv_min= sigma0_inv(M_min);
  sinv_max= sigma0_inv(M_max);
}

Sigma::~Sigma()
{
  gsl_interp_free(interp_);
  gsl_interp_free(interp2_);
  gsl_interp_accel_free(acc_);
  gsl_interp_accel_free(acc2_);

  free(M_);
}

double Sigma::M(const double sigma0) const
{
  // sigma0 is sigma(M, z=0)
  // return value: M(sigma)
#ifdef DEBUG
  if(1/sigma0 < sinv_[0] || 1/sigma0 > sinv_[n-1]) {
    cerr << "Error: sigma_M(sigma0) out of range.\n";
    cerr << sinv_[0] << " " << sinv_[n-1] << " " << 1.0/sigma0 << endl;
  }
#endif
  return gsl_interp_eval(interp_, sinv_, M_, 1.0/sigma0, acc_);
}

double Sigma::sigma0_inv(const double MM) const
{
  return gsl_interp_eval(interp2_, M_, sinv_, MM, acc2_);
}

