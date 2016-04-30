#include <iostream>
#include <cassert>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "msg.h"
#include "hod.h"
#include "catalogue.h"
#include "satellite.h"

#include "sky.h"
#include "mf_cumulative.h"
using namespace std;


static gsl_rng* rng= 0;


//
// Class member functions
//
Catalogue::Catalogue() :
  tree(0), ntree(0), ncen(0), nsat(0)
{
  
}

Catalogues::Catalogues()
{

}

Catalogues::Catalogues(const size_t n)
{
  allocate(n);
}

Catalogues::~Catalogues()
{
  for(Catalogues::iterator p= begin(); p != end(); ++p) {
    delete *p;
  }
}

void Catalogues::allocate(const size_t n)
{
  assert(empty());
  
  for(int i=0; i<n; ++i)
    push_back(new Catalogue());

  msg_printf(msg_verbose, "allocated %d new catalogues\n");

  assert(size() == n);
}


//
// catalogue module
//
void catalogue_init()
{
  if(rng)
    return;
  
  const unsigned int seed= (unsigned int) time(NULL);
  rng= gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(rng, seed);

  satellite_init(rng);

  msg_printf(msg_info, "# catalogue random initialised with %u\n", seed);
}

void catalogue_free()
{
  satellite_free();
  gsl_rng_free(rng);
}


//
// Generate a mock catalogue from a lightcone
//


void catalogue_generate_mock(Hod* const hod,
			     LightCone const * const lightcone,
			     const double z_min, const double z_max,
			     Catalogue * const cat)
{
  assert(hod);
  assert(lightcone);
  assert(cat);

  cat->clear();
  
  Particle p; p.w= 1.0;
  int ncen_total= 0, nsat_total= 0;


  cerr << "src lightcone " << lightcone->size() << endl;
  cerr << "target catalogue" << cat->size() << endl;
  
  for(LightCone::const_iterator h=
	lightcone->begin(); h != lightcone->end(); ++h) {
    if(h->z < z_min || h->z >= z_max)
      continue;

    hod->compute_param_z(h->z);

    // centrals
    double ncen= hod->ncen(h->M);

    if(gsl_rng_uniform(rng) > ncen)
      continue;

    p.x[0]= h->x[0];
    p.x[1]= h->x[1];
    p.x[2]= h->x[2];
    p.vr  = h->v[0];
    p.z   = h->z;
    p.radec[0] = h->radec[0];
    p.radec[1] = h->radec[1];

    cat->push_back(p);
    ncen_total++;

    double nsat_mean= hod->nsat(h->M);
    int nsat= gsl_ran_poisson(rng, nsat_mean);
    nsat_total += nsat;

    // satellites
    for(int isat=0; isat<nsat; ++isat) {
      satellite(&*h, &p);
      p.x[0] += h->x[0];
      p.x[1] += h->x[1];
      p.x[2] += h->x[2];
      p.z     = h->z;
      p.vr   += h->v[0];
      p.radec[0] = h->radec[1];
      p.radec[1] = h->radec[1];

      cat->push_back(p);
    }    
  }
  cat->ncen= ncen_total;
  cat->nsat= nsat_total;
}

/*
void catalogue_generate_randoms(Hod* const hod,
				Sky* const sky,
				MfCumulative* const mc,
				Catalogue * const cat)
{
  Particle particle;
  Particle* p= &particle;

  const double vol= sky->width[0]*sky->width[1]*sky->width[2];
  const int nhalo= vol*mc->nM_max;

  // Centrals
  for(int i=0; i<nhalo; ++i) {
    for(int k=0; k<3; ++k)
      p->x[k]= sky->left[k] + sky->width[k]*gsl_rng_uniform(rng);

    sky_compute_ra_dec(sky, p);
    if(p->ra < sky->ra_range[0] ||  p->ra > sky->ra_range[1] ||
       p->dec < sky->dec_range[0] || p->dec > sky->dec_range[1])
      continue;

    double nM= mc->nM_max*gsl_rng_uniform(rng);
    double M= mc->M(nM);

    double ncen= hod->ncen(h->M);
    //if(gsl_rng_uniform(rng) > ncen)

  }

  
}
*/
