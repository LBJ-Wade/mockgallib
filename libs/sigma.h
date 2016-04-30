#ifndef SIGMA_H
#define SIGMA_H 1

#include <gsl/gsl_interp.h>
#include "power.h"

class Sigma {
 public:
  Sigma(PowerSpectrum const * const ps,
	const double M_min=1.0e10, const double M_max=1.0e16, const int n=1001);
  ~Sigma();
  double sigma0_inv(const double M);
  double M(const double sigma0);

  int n;
  double M_min, M_max;
  double sinv_min, sinv_max;
  double *M_, *sinv_;
 private:
  gsl_interp *interp_, *interp2_;
  gsl_interp_accel *acc_, *acc2_;
};

#endif
