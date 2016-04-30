#include "util.h"
#include "sky.h"
#include "lightcone.h"
#include "remap.h"
#include "slice.h"
#include "cola_file.h"
#include "distance.h"
#include "halo_mass.h"
#include "halo_concentration.h"



LightCones::~LightCones()
{
  for(LightCones::iterator p= begin(); p != end(); ++p) {
    delete *p;
  }
}

void fill_lightcone(const char filename[],
		    const float a_snp,
		    const float r_min, const float r_max,
		    Sky const * const sky,
		    Remap const * const remap,
		    Slice const * const slice,
		    HaloMassFoF const * const halo_mass,
		    LightCone* const lightcone)
{
  // halo_concentration_init is required before using this function
  
  cola_halo_file_open(filename);

  Halo halo;
  Halo* const h= &halo;

  while(cola_halo_file_read_one(h)) {
    // remap cube to cuboid
    if(!remap->coordinate(h))
      continue;

    // cuboid to sliced cuboid (multiple mock from one cuboid)
    h->slice= slice->transform(h->x);

    h->a= a_snp;
    h->r= util::norm(h->x);

    if(h->r < r_min || h->r >= r_max)
	continue;

    // compute redshift at halo position
    h->z= distance_redshift(h->r);
      
    // compute sky position
    sky->compute_radec(h->x, h->radec);

    if(h->radec[0] < sky->ra_range[0]  || h->radec[0] > sky->ra_range[1] ||
       h->radec[1] < sky->dec_range[0] || h->radec[1] > sky->dec_range[1])
      continue;

    // convert nfof to halo mass
    h->M= halo_mass->mass(h->nfof);

    // set halo concentration / rs
    h->rs= halo_concentration_rs(h);
      
    lightcone->push_back(*h);
  }
   
  cola_halo_file_close();
}
