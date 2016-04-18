#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include "cosmology.h"

static double D0= 0.0;  // normalization factor

static double growth_int(double a, void* param);
static double growth(const double a);

double growth_D(const double a)
{
  // Growth factor D(a)
  if(D0 == 0.0)
    D0= growth(1.0);

  if(a == 0.0)
    return 0.0;

  return growth(a)/D0;
}

double growth_f(const double a)
{
  const double om= cosmology_omega_m();
  
  if(D0 == 0.0)
    D0= growth(1.0);    
  
  // Growth rate f=dlnD/dlna
  double hubble_a= sqrt(om/(a*a*a) + 1.0 - om);
  double d= growth(a)/D0;

  double f_ex= 1.0/(d*D0*a*a*hubble_a*hubble_a)
                 - 1.5*om/(hubble_a*hubble_a*a*a*a);

  return f_ex;
}



double growth_integrand(double a, void* params)
{
  const double om= *(double*) params;
  return pow(a/(om + (1.0 - om)*a*a*a), 1.5);
}

double growth(const double a)
{
  // Compute integral of D_tilda(a)
  // growth_factor = D_tilda(a)/D_tilda(1.0)

  const double om= cosmology_omega_m();
  const double hubble_a= sqrt(om/(a*a*a) + 1.0 - om);

  const int worksize= 1000;
  double result, abserr;

  gsl_integration_workspace *workspace=
    gsl_integration_workspace_alloc(worksize);

  gsl_function F;
  F.function = &growth_integrand;
  F.params= (void*) &om;

  gsl_integration_qag(&F, 0, a, 0, 1.0e-8, worksize, GSL_INTEG_GAUSS41,
		      workspace, &result, &abserr);

  gsl_integration_workspace_free(workspace);

  return hubble_a * result;
}
