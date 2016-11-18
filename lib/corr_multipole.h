#ifndef CORR_MULTIPOLE_H
#define CORR_MULTIPOLE_H 1

#include "hist2d.h"
#include "catalogue.h"

struct CorrMultipole {
  CorrMultipole(const int nbin);
  ~CorrMultipole();
  int n;
  double* r;
  double* xi0;
  double* xi2;
};

void corr_multipole_init(const float r_min_, const float r_max_, const int nbin_, const int nbin_mu_);
void corr_multipole_free();
void corr_multipole_set_pair_correction(const char filename[]);
void corr_multipole_set_radec_min(const float ra_min_, const float dec_min_);
CorrMultipole* corr_multipole_i(const int i);
void corr_multipole_compute_pairs_rr(Catalogues* const cats_rand,
				     Histogram2D<LogBin, LinearBin>* const rr);
void corr_multipole_compute_with_rr(Catalogues* const cats_data,
				    Catalogues* const cats_rand,
				    Histogram2D<LogBin, LinearBin> const * const rr);
void corr_multipole_compute(Catalogues* const cats_data,
			    Catalogues* const cats_rand);


#endif
