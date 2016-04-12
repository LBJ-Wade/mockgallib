#ifndef TINKER_H
#define TINKER_H 1

#include <gsl/gsl_integration.h>

struct MF {
  double z;
  double alpha;
  gsl_integration_cquad_workspace* w;
};

MF* mf_alloc();
void mf_free(MF* const mf);
void mf_set_redshift(MF* const mf, const double a);
double mf_f(MF const * const mf, const double nu);
double mf_b(const double nu);

#endif

