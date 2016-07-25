#ifndef SIGMA_H
#define SIGMA_H 1

#include <gsl/gsl_interp.h>
#include "power.h"


void sigma_init(const double M_min=1.0e10, const double M_max=1.0e16,
	   const int n=1001);
void sigma_free();
double sigma_M(const double sigma0);
double sigma_inv(const double M);

double* sigma_M_array();
double* sigma_sinv_array();
double sigma_M_min();
double sigma_M_max();
double sigma_sinv_min();
double sigma_sinv_max();
int sigma_n();
#endif
