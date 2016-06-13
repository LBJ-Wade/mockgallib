#ifndef CORR_PROJECTED_H
#define CORR_PROJECTED_H 1

#include <vector>
#include "particle.h"
#include "kdtree.h"

struct CorrProjected {
  CorrProjected(const int nbin);
  ~CorrProjected();
  int n;
  double* rp;
  double* wp;
  double* dwp;
};

/*
struct Catalogue: std::vector<Particle>
{
  KDTree* tree;
  size_t ntree;
  int ncen, nsat;
};


void corr_projected_init(const float rp_min_, const float rp_max_, const int nbin_, const float pi_max_, const int pi_nbin_);
void corr_projected_free();

void corr_projected_compute(std::vector<Catalogue>& vdata,
			    std::vector<Catalogue>& vrandom,
			    std::vector<CorrProjected*>& vcorr);

void corr_alloc(const int n_data_cat, const int nbin, std::vector<CorrProjected*>& vcorr);

void corr_projected_write(const int index, const std::vector<CorrProjected*>& vcorr);

void corr_projected_summarise(const std::vector<CorrProjected*>& vcorr,
			      CorrProjected* corr_out);

*/
  

#endif
