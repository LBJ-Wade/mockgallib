#ifndef HALO_MASS_H
#define HALO_MASS_H 1

#include <gsl/gsl_spline.h>

class HaloMassFoF {
 public:
  HaloMassFoF(const char filename[], const double M0=1.0e14);
  ~HaloMassFoF();
  double mass(const int nfof) {
    // statistics if poor for large rare haloes
    if(nfof > nfof0_) return m_ * nfof;
  
    return gsl_spline_eval(spline_, nfof, acc_);
  }
 private:
  gsl_spline* spline_;
  gsl_interp_accel* acc_;
  double nfof0_, M0_, m_;
};

#endif
