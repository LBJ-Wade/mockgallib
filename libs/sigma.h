#ifndef SIGMA_H
#define SIGMA_H 1

#include <gsl/gsl_interp.h>
#include "power.h"

struct Sigma {
  int n; //= 1001;
  double *M, *sinv;
  double M_min, M_max;
  gsl_interp *interp, *interp2;
  gsl_interp_accel *acc, *acc2;
  double sinv_min, sinv_max;
};

Sigma* sigma_alloc(PowerSpectrum const * const ps,
		   const double M_min=1.0e10, const double M_max=1.0e16, const int n=1001);
void sigma_free(Sigma* const s);
double sigma_M(Sigma* const s, const double sigma0);
double sigma_0inv(Sigma* const s, const double M);
  
#endif
