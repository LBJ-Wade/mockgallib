#ifndef COSMOLOGY_H
#define COSMOLOGY_H

void cosmology_set_omega_m(const double omega_m_);

double cosmology_omega_m();
double cosmology_rho_m();

double cosmology_compute_comoving_distance(const double a);

static inline double cosmology_redshift(const double a)
{
  return 1.0/a - 1.0;
}

static inline double cosmology_scalefactor(const double z)
{
  return 1.0/(1.0 + z);
}

#endif
