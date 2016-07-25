#include <cassert>
#include <time.h>
#include "rand.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

static gsl_rng* rng= 0;

void rand_init()
{
  if(rng) return;
  
  const unsigned int seed= (unsigned int) time(NULL);
  rng= gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(rng, seed);
}


void rand_free()
{
  if(rng)
    gsl_rng_free(rng);
}


double rand_uniform()
{
#ifdef DEBUG
  assert(rng);
#endif
  return gsl_rng_uniform(rng);
}


double rand_poisson(const double lmbda)
{
#ifdef DEBUG
  assert(rng);
#endif

  return gsl_ran_poisson(rng, lmbda);
}

double rand_gaussian()
{
#ifdef DEBUG
  assert(rng);
#endif

  return gsl_ran_ugaussian(rng);
}
