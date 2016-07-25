#ifndef NFW_H
#define NFW_H 1

#include "power.h"
#include "sigma.h"
#include "halo.h"

void halo_concentration_init();

float halo_concentration_rs(Halo const * const h);

#endif
