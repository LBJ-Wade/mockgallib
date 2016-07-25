#ifndef SATELLITE_H
#define SATELLITE_H 1

#include <gsl/gsl_rng.h>
#include "particle.h"


void satellite_init();
void satellite_free();
void satellite(Halo const * const h, Particle* const g);

#endif
