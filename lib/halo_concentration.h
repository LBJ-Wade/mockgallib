#ifndef NFW_H
#define NFW_H 1

#include "power.h"
#include "sigma.h"
#include "halo.h"

void halo_concentration_init(PowerSpectrum const * const ps);
void halo_concentration_init(Sigma const * const s);

float halo_concentration_rs(Halo const * const h);

#endif
