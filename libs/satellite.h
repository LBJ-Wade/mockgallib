#ifndef SATELLITE_H
#define SATELLITE_H 1

#include <gsl/gsl_rng.h>
#include "particle.h"


void satellite_init(gsl_rng* const rng);
void satellite_free();
void satellite(Halo const * const h, Particle* const g);

#endif
