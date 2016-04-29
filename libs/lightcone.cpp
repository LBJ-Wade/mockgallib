#include "util.h"
#include "sky.h"
#include "lightcone.h"
#include "remap.h"
#include "slice.h"
#include "cola_file.h"
#include "distance.h"

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
		    LightCone* const lightcone)
{
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
    sky_compute_ra_dec(sky, h);

    if(h->ra < sky->ra_range[0] || h->ra > sky->ra_range[1] ||
       h->dec < sky->dec_range[0] || h->dec > sky->dec_range[1])
      continue;

    // convert nfof to halo mass
    // ToDo!!!
    //h->M= halo_mass(h->nfof);

    // set halo concentration / rs
    // ToDo !!!
    //h->rs= compute_halo_rs(h);
      
    lightcone->push_back(*h);
  }
   
  cola_halo_file_close();
}
