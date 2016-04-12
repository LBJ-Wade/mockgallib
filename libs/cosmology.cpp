#include "cosmology.h"
#include "const.h"

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

