#ifndef NBAR_H
#define NBAR_H 1

#include <vector>

#include "power.h"
#include "sigma.h"
#include "mf.h"
#include "hod.h"

struct Nbar {
  double z, nbar, dnbar, ncen, nsat;
};

struct NbarIntegration {
  NbarIntegration();
  NbarIntegration(NbarIntegration const * const ni);
  gsl_integration_cquad_workspace* w;
  MF* mf;
  Hod* hod;
  double rho_m;
  double D;
  double z;
};


NbarIntegration* nbar_integration_alloc(Hod* const hod);
void nbar_integration_free(NbarIntegration* const ni);
double nbar_compute(NbarIntegration* const ni, const double z);

#endif
