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
  gsl_integration_cquad_workspace* w;
  Sigma* s;
  MF* mf;
  Hod* hod;
  double rho_m;
  double D;
  double z;
};


NbarIntegration* nbar_integration_alloc(PowerSpectrum const * const ps,
					Hod* const hod);
void nbar_integration_free(NbarIntegration* const ni);
double nbar_compute(NbarIntegration* const ni, const double z);


//void nbar_init(PowerSpectrum const * const ps, const double omega_m, const double z_max);
//void nbar_free();

//void nbar_read(const char filename[], const double z_min, const double z_max,
//	       std::vector<Nbar>& v);

//void nbar_compute_nz(const double c[], std::vector<Nbar>* const v);
//void nbar_compute_derivatives_dnz(const double c[], HodParam* const hod,
//				  std::vector<Nbar>& v);

//void nbar_compute(HodParam* const hod, std::vector<Nbar>* const v);

//double compute_n_hod_derivative(HodParam const * const hod, const int j);
//void nbar_compute_derivatives(HodParam* const hod, std::vector<Nbar>& v);
#endif
