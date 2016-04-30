#include <cstdio>
#include <cmath>
#include <cassert>
#include <gsl/gsl_rng.h>

#include "msg.h"
#include "cosmology.h"
#include "halo_concentration.h"


static double M_star= 0.0;

void halo_concentration_init(PowerSpectrum const * const ps)
{
  // Compute Mstar for cmean
  Sigma s(ps);
  M_star= s.M(1.0);
  msg_printf(msg_verbose, "Mstar= %e\n", M_star);
}

void halo_concentration_init(Sigma const * const s)
{
  M_star= s->M(1.0);
  msg_printf(msg_verbose, "Mstar= %e\n", M_star);
}

static inline double cmean(const double M, const double a)
{
  // Bullock concentraion
  //return 11.0*a*pow(M200/2.344228815319e12, -0.13);
#ifdef DEBUG
  assert(M_star > 0.0);
#endif

  return 11.0*a*pow(M/M_star, -0.13);

  // 2.34e12 /h M_solar is M* at z=0 which depends on cosmology
  // ToDo: check M* is correct
}

float halo_concentration_rs(Halo const * const h)
{
  double a= h->a;
  double rho_m= cosmology_rho_m()/(a*a*a);
  double r200m= pow(h->M/(4.0*M_PI/3.0*200.0*rho_m), 1.0/3.0);
    // physical 1/h Mpc
  
  return r200m / cmean(h->M, h->a);
}

