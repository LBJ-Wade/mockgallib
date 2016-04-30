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

Sigma* sigma_alloc(PowerSpectrum const * const ps,
		   const double M_min, const double M_max, const int n)
{
  const double rho_m= cosmology_rho_m(); assert(rho_m >= 0.0);

  Sigma* s= new Sigma();
  s->M_min= M_min;
  s->M_max= M_max;
  s->n= n;
  const double log_Mmin= 0.99*log(M_min);
  const double log_Mmax= 1.01*log(M_max);

  s->M= (double*) malloc(sizeof(double)*2*n); assert(s->M);
  s->sinv= s->M + n;
  
  for(int i=0; i<n; ++i) {
    double logM = log_Mmin + (log_Mmax - log_Mmin)*((double) i/(n-1));
    s->M[i]= exp(logM);
    double R= pow(s->M[i]/(4.0*M_PI/3.0*rho_m), 1.0/3.0);

    s->sinv[i]= 1.0/ps->compute_sigma(R);
  }

  // Function: 1/sigma -> M
  s->interp= gsl_interp_alloc(gsl_interp_cspline, n);
  gsl_interp_init(s->interp, s->sinv, s->M, n);
  s->acc= gsl_interp_accel_alloc(); 

  // Function: M -> 1/sigma
  s->interp2= gsl_interp_alloc(gsl_interp_cspline, n);
  gsl_interp_init(s->interp2, s->M, s->sinv, n); 
  s->acc2= gsl_interp_accel_alloc();

  s->sinv_min= sigma_0inv(s, M_min);
  s->sinv_max= sigma_0inv(s, M_max);

  return s;
}

void sigma_free(Sigma* const s)
{
  gsl_interp_free(s->interp);
  gsl_interp_free(s->interp2);
  gsl_interp_accel_free(s->acc);
  gsl_interp_accel_free(s->acc2);

  free(s->M);

  delete s;
}

double sigma_M(Sigma* const s, const double sigma0)
{
  // sigma0 is sigma(M, z=0)
  // return value: M(sigma)
#ifdef DEBUG
  if(1/sigma0 < s->sinv[0] || 1/sigma0 > s->sinv[s->n-1]) {
    cerr << "Error: sigma_M(sigma0) out of range.\n";
    cerr << s->sinv[0] << " " << s->sinv[s->n-1] << " " << 1.0/sigma0 << endl;
  }
#endif
  return gsl_interp_eval(s->interp, s->sinv, s->M, 1.0/sigma0, s->acc);
}

double sigma_0inv(Sigma* const s, const double M)
{
  return gsl_interp_eval(s->interp2, s->M, s->sinv, M, s->acc2);
}

/*
double sigma_0inv_min()
{
  return sinv_min; //[0];
}

double sigma_0inv_max()
{
  return sinv_max;//return sinv[n-1];
}
*/
