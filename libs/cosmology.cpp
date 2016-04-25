#include "cosmology.h"
#include "const.h"

#include <gsl/gsl_integration.h>

static double omega_m= 0.3;  // Omega_m(z=0)

void cosmology_set_omega_m(const double omega_m_)
{
  omega_m= omega_m_;
}

double cosmology_omega_m()
{
  return omega_m;
}

double cosmology_rho_m()
{
  // Mean matter density at z=0, or comoving mean matter density
  return omega_m*rho_crit_0;
}

static double distance_integrand(double a, void* params) {
  // 1/[ a^2 H(a)/H0 ]
  const double om= *(double*) params;
  return 1.0/sqrt(om*a + (1.0 - om)*(a*a*a*a));
}

double cosmology_compute_comoving_distance(const double a)
{
  // Light emitted at comoving distance r at scale facgor a reach
  // the observer now
  //
  // r = \int c*dt/aH(a) = \int c*da/[a^2 H(a)]
  //   = c/H0 \int da/[a^2 sqrt(omega_m*a + omega_lambda*a^4)]
  
  const int nwork= 1000;
  gsl_integration_workspace* w= gsl_integration_workspace_alloc(nwork);

  gsl_function F;
  F.function= &distance_integrand;
  F.params= (void*) &omega_m;

  double result, error;
  gsl_integration_qags(&F, a, 1.0, 0, 1.0e-8, nwork, w, &result, &error);

  gsl_integration_workspace_free(w);

  return cH0inv*result;
}
