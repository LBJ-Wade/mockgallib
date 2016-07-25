//
// Generate a lightcone from snapshots
//

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
#include "rand.h"

//
// Dependence:
//     sigma_init() - power_init()
//     distance(z_max)
//
using namespace std;

static void fill_lightcone(Snapshot const * const snp,
			   Sky const * const sky,
			   Remap const * const remap,
			   Slice const * const slice,
			   const bool random,
			   LightCones* const lightcones);

static inline void
  randomise_position(const float boxsize, float* const x)
{
  x[0]= boxsize*rand_uniform();
  x[1]= boxsize*rand_uniform();
  x[1]= boxsize*rand_uniform();
}

void cola_lightcones_create(Snapshots const * const snapshots,
			    Sky const * const sky,
			    Remap const * const remap,
			    Slice const * const slice,
			    LightCones* const lightcones,
			    const bool random)
{
  // Prerequisite: sigma_init()
  halo_concentration_init();
  distance_init(sky->z_range[1]);
      
  if(random) {
    rand_init();
    msg_printf(msg_verbose, "Generate random lightcone\n");
  }
  else
    msg_printf(msg_verbose, "Generate halo lightcone\n");    
  
  lightcones->clear();
  
  for(Snapshots::const_iterator snp= snapshots->begin();
      snp != snapshots->end(); ++snp) {
    
    fill_lightcone(*snp, sky, remap, slice, random, lightcones);
  }
}	   

void fill_lightcone(Snapshot const * const snp,
		    Sky const * const sky,
		    Remap const * const remap,
		    Slice const * const slice,
		    const bool random,
		    LightCones* const lightcones)
{
  // Prerequisite:
  //   
  // halo_concentration_init() -- sigma_init() -- power_init()
  //

  if(lightcones->size() < slice->n) {
    lightcones->resize(slice->n);
  }

  msg_printf(msg_verbose, "filling lightcone from %s, a=%.3f\n",
	     snp->filename, snp->a_snp);

  float boxsize;
  
  cola_halo_file_open(snp->filename, &boxsize);
  // May throw ColaFileError()

  Halo halo;
  Halo* const h= &halo;

  const float r_min= snp->r_min;
  const float r_max= snp->r_max;
  const float ra_min= sky->ra_range[0];
  const float ra_max= sky->ra_range[1];
  const float dec_min= sky->dec_range[0];
  const float dec_max= sky->dec_range[1];

  while(cola_halo_file_read_one(h)) {
    // randomise the coordinate if random = true
    if(random)
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
   
  cola_halo_file_close();
}
