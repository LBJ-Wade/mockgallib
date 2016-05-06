#include <vector>

#include "util.h"
#include "sky.h"
#include "lightcone.h"
#include "remap.h"
#include "slice.h"
#include "cola_file.h"
#include "distance.h"
#include "halo_mass.h"
#include "halo_concentration.h"
#include "snapshot.h"

using namespace std;

static void fill_lightcone(Snapshot const * const snp,
			   Sky const * const sky,
			   Remap const * const remap,
			   Slice const * const slice,
			   LightCones* const lightcones);


void cola_lightcone_create(Snapshots const * const snapshots,
			   Sky const * const sky,
			   Remap const * const remap,
			   Slice const * const slice,
			   LightCones* const lightcones)
{
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
  
  cola_halo_file_open(snp->filename);

  Halo halo;
  Halo* const h= &halo;

  const float r_min= snp->r_min;
  const float r_max= snp->r_max;
  const float ra_min= sky->ra_range[0];
  const float ra_max= sky->ra_range[1];
  const float dec_min= sky->dec_range[0];
  const float dec_max= sky->dec_range[1];

  while(cola_halo_file_read_one(h)) {
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
