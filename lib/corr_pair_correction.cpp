#include <iostream>
#include <cstdio>
#include <cassert>
#include <gsl/gsl_spline.h>

#include "corr_pair_correction.h"

using namespace std;

static gsl_spline* spline= 0;
static gsl_interp_accel* acc= 0;
static double theta_max= 0;
static double theta_min= 0;

void corr_pair_correction_init(const char filename[])
{
  const int n_alloc= 50;
  
  double* theta= (double*) malloc(sizeof(double)*n_alloc*2); assert(theta);
  double* w= theta + n_alloc;

  FILE* fp= fopen(filename, "r");

  char buf[128];
  int n= 0;
  while(fgets(buf, 127, fp)) {
    assert(n < n_alloc);
    sscanf(buf, "%le %le", theta+n, w+n);
    n++;
  }

  spline= gsl_spline_alloc(gsl_interp_cspline, n);
  acc= gsl_interp_accel_alloc();

  gsl_spline_init(spline, theta, w, n);
  theta_min= theta[0];
  theta_max= theta[n-1];

  fclose(fp);
  free(theta);
}

void corr_pair_correction_free()
{
  if(spline)
    gsl_spline_free(spline);
  if(acc)
    gsl_interp_accel_free(acc);
}

double corr_pair_correction(const double theta)
{
  if(theta >= theta_max)
    return 1.0f;

  if(theta <= theta_min) {
    //cerr << "theta_min " << theta << " < " << theta_min << endl;
    return -1.0;
  }

  assert(spline);
  return gsl_spline_eval(spline, theta, acc);
}
