#ifndef CORR_PROJECTED_H
#define CORR_PROJECTED_H 1

#include <vector>
#include "particle.h"
#include "catalogue.h"
#include "hist2d.h"

struct CorrProjected {
  CorrProjected(const int nbin);
  ~CorrProjected();
  void print(FILE* fp);
  int n;
  double* rp;
  double* wp;
  double* dwp;
};

void corr_projected_init(const float rp_min_, const float rp_max_, const int nbin_, const float pi_max_, const int nbin_pi_);

void set_radec_min(const float ra_min_, const float dec_min_);
void corr_set_pair_correction(const char filename[]);

void corr_projected_compute(Catalogues* const cats_data,
			    Catalogues* const cats_rand,
			    CorrProjected* const corr);

void corr_projected_compute_pairs_rr(Catalogues* const cats_rand,
				     Histogram2D<LogBin, LinearBin>* const rr);
void corr_projected_compute_pairs_all(Catalogues* const cats_data,
				      Catalogues* const cats_rand,
				      Histogram2D<LogBin, LinearBin>* const dd,
				      Histogram2D<LogBin, LinearBin>* const dr,
				      Histogram2D<LogBin, LinearBin>* const rr);
void corr_projected_compute_with_rr(Catalogues* const cats_data,
				    Catalogues* const cats_rand,
			    Histogram2D<LogBin, LinearBin> const * const rr);

CorrProjected* corr_projected_i(const int i);

void corr_projected_compute_direct(Catalogues* const cats_data,
				   Catalogues* const cats_rand,
				   CorrProjected* const corr);
void corr_projected_compute_pairs_rr_direct(Catalogues* const cats_rand,
				   Histogram2D<LogBin, LinearBin>* const rr);
void corr_projected_compute_pairs_all_direct(Catalogues* const cats_data,
					     Catalogues* const cats_rand,
				     Histogram2D<LogBin, LinearBin>* const dd,
				     Histogram2D<LogBin, LinearBin>* const dr,
				     Histogram2D<LogBin, LinearBin>* const rr);

/*
struct Catalogue: std::vector<Particle>
{
  KDTree* tree;
  size_t ntree;
  int ncen, nsat;
};


void corr_projected_init(const float rp_min_, const float rp_max_, const int nbin_, const float pi_max_, const int pi_nbin_);
void corr_projected_free();


void corr_alloc(const int n_data_cat, const int nbin, std::vector<CorrProjected*>& vcorr);

void corr_projected_write(const int index, const std::vector<CorrProjected*>& vcorr);

void corr_projected_summarise(const std::vector<CorrProjected*>& vcorr,
			      CorrProjected* corr_out);

*/
  

#endif
