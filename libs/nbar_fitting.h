#ifndef NBAR_FITTING_H
#define NBAR_FITTING_H 1

#include "nbar.h"

struct NbarFitting {
  Hod* hod;
  std::vector<NbarIntegration*> vni;
  double z_min, z_max;
  std::vector<Nbar> const * vobs;
  std::vector<Nbar>* vhod;
};

NbarFitting* nbar_fitting_alloc(PowerSpectrum const * const ps,
				Hod* const hod,
				const std::vector<Nbar>* const vnbar_obs,
				const double z_min, const double z_max);
void nbar_fitting_free(NbarFitting* const fitting);
void nbar_fitting_compute(NbarFitting* fitting);

#endif
