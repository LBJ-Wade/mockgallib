#include <iostream>
#include <vector>
#include <queue>
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
// Satelite queue class
//

struct SatelliteRandom {
  double M;
  int n;
};
  
class SatelitteQueue {
 public:
  bool empty() const {
    return q.empty();
  }
  double pop() {
    double M= q.front().M;
    if(q.front().n == 1)
      q.pop();
    else
      q.front().n--;
    return M;
  }
  void push(const double M, const int n) {
    SatelliteRandom s;
    s.M= M; s.n= n;
    q.push(s);
  }
      
 private:
  queue<SatelliteRandom> q;
};


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

  cerr << ncen_total << " centrals, " << nsat_total << " statellites\n";
}


void catalogue_generate_centrals(Hod* const hod,
      LightCone const * const lightcone, const double z_min, const double z_max,
      Catalogue * const cat)
{
  assert(hod);
  assert(lightcone);
  assert(cat);

  cat->clear();
  
  Particle p; p.w= 1.0;
  int ncen_total= 0, nsat_total= 0;

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

  }
  cat->ncen= ncen_total;
  cat->nsat= 0;
}


void catalogue_generate_random(Hod* const hod,
			       LightCone const * const lightcone,
			       const double z_min, const double z_max,
			       Catalogue * const cat)
{
  // lightcone: random lightcone
  assert(hod);
  assert(lightcone);
  assert(cat);

  cat->clear();
  
  Particle p; p.w= 1.0;
  int ncen_total= 0, nsat_total= 0;

  const double dz= 0.02;
  const int n_zbin= (int) ceil((z_max - z_min)/dz);
  vector<SatelitteQueue> qsat(n_zbin);
  
  for(LightCone::const_iterator h=
	lightcone->begin(); h != lightcone->end(); ++h) {
    if(h->z < z_min || h->z >= z_max)
      continue;

    hod->compute_param_z(h->z);

    // centrals
    double ncen= hod->ncen(h->M);
    int iz= (int)((h->z - z_min)/dz); assert(0 <= iz && iz < n_zbin);
    
    if(gsl_rng_uniform(rng) <= ncen) {
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

      qsat[iz].push(h->M, nsat);
      nsat_total += nsat;
    }
    else if(!qsat[iz].empty()) {
      // put satellites
      Halo hh= *h;
      hh.M= qsat[iz].pop();
      satellite(&hh, &p);
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
