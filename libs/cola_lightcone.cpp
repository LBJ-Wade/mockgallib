#include <cstdio>
#include <iostream>
#include <vector>

#include "msg.h"
#include "util.h"
#include "cola_file.h"
#include "distance.h"
#include "halo_mass.h"
#include "halo_concentration.h"
#include "cola_lightcone.h"

using namespace std;
static gsl_rng* rng= 0;

static void fill_lightcone(Snapshot const * const snp,
			   Sky const * const sky,
			   Remap const * const remap,
			   Slice const * const slice,
			   LightCones* const lightcones);

static inline void
  randomise_position(const float boxsize, float* const x)
{
  x[0]= boxsize*gsl_rng_uniform(rng);
  x[1]= boxsize*gsl_rng_uniform(rng);
  x[1]= boxsize*gsl_rng_uniform(rng);
}

void cola_lightcones_create(Snapshots const * const snapshots,
			    Sky const * const sky,
			    Remap const * const remap,
			    Slice const * const slice,
			    LightCones* const lightcones,
			    gsl_rng* const random_generator)
{
  rng= random_generator;
  lightcones->clear();
  
  for(Snapshots::const_iterator snp= snapshots->begin();
      snp != snapshots->end(); ++snp) {
    
    fill_lightcone(*snp, sky, remap, slice, lightcones);
  }
}	   

void fill_lightcone(Snapshot const * const snp,
		    Sky const * const sky,
		    Remap const * const remap,
		    Slice const * const slice,
		    LightCones* const lightcones)
{
  // halo_concentration_init is required before using this function

  if(lightcones->size() < slice->n) {
    lightcones->resize(slice->n);
  }

  msg_printf(msg_verbose, "filling lightcone from %s, a=%.3f\n",
	     snp->filename, snp->a_snp);

  float boxsize;
  
  cola_halo_file_open(snp->filename, &boxsize);

  Halo halo;
  Halo* const h= &halo;

  const float r_min= snp->r_min;
  const float r_max= snp->r_max;
  const float ra_min= sky->ra_range[0];
  const float ra_max= sky->ra_range[1];
  const float dec_min= sky->dec_range[0];
  const float dec_max= sky->dec_range[1];

  while(cola_halo_file_read_one(h)) {
    // randomise the coordinate if rng != 0
    if(rng)
      randomise_position(boxsize, h->x);

    // remap cube to cuboid
    if(!remap->coordinate(h))
      continue;

    // cuboid to sliced cuboid (multiple mock from one cuboid)
    slice->transform(h);

    h->a= snp->a_snp;
    h->r= util::norm(h->x);

    if(h->r < r_min || h->r >= r_max)
	continue;

    // compute sky ra-dec
    sky->compute_radec(h->x, h->radec);

    if(h->radec[0] < ra_min  || h->radec[0] > ra_max ||
       h->radec[1] < dec_min || h->radec[1] > dec_max)
      continue;

    // compute raidial velocity
    h->vr= util::dot(h->x, h->v)/util::norm(h->x);

    // compute redshift at halo position
    h->z= distance_redshift(h->r);

    // convert nfof to halo mass
    h->M= snp->halo_mass->mass(h->nfof);

    // set halo concentration / rs
    h->rs= halo_concentration_rs(h);

    (*lightcones)[h->slice]->push_back(*h);
  }

  /* ???
  for(LightCones::iterator p=
	lightcones->begin(); p != lightcones->end(); ++p) {
  }
  */
   
  cola_halo_file_close();
}
