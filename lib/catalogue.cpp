#include <iostream>
#include <vector>
#include <queue>
#include <cassert>

#include "msg.h"
#include "hod.h"
#include "rand.h"
#include "util.h"
#include "catalogue.h"
#include "satellite.h"

#include "sky.h"
#include "mf_cumulative.h"
using namespace std;

static bool initialised= false;

//
// Satelite queue class
//

struct SatelliteRandom {
  float M, rs;
  int n;
};
  
class SatelitteQueue {
 public:
  bool empty() const {
    return q.empty();
  }
  void pop(Halo& h) {
    assert(q.front().n > 0);
    h.M= q.front().M;
    h.rs= q.front().rs;


    if(q.front().n == 1) {
      q.pop();
    }
    else {
      q.front().n--;
    }
  }
  void push(const float M, const float rs, const int n) {
    SatelliteRandom s;
    s.M= M; s.rs= rs; s.n= n;
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
  if(!initialised) {
    rand_init();
    satellite_init();

    msg_printf(msg_debug, "catalogue module initialised\n");
    initialised= true;    
  }
}

void catalogue_free()
{
  if(initialised)
    satellite_free();
}


//
// Generate a mock catalogue from a lightcone
//

void catalogue_generate_mock(Hod* const hod,
			     LightCone const * const lightcone,
			     const double z_min, const double z_max,
			     Catalogue * const cat)
{
  catalogue_init();
  assert(hod);
  assert(lightcone);
  assert(cat);

  cat->clear();
  
  Particle p; p.w= 1.0;
  int ncen_total= 0, nsat_total= 0;

  /*
  cerr << "src lightcone " << lightcone->size() << endl;
  cerr << z_min << " " << z_max << endl;
  cerr << "HOD c " << hod->c[0] << " " << hod->c[1] << " " << hod->c[2] << " "
       << hod->c[3] << " " << hod->c[4] << " " << hod->c[5] << " "
       << hod->c[6] << " " << hod->c[7] << " " << hod->c[8] << " " << hod->c[9]
       << endl;
  */


  for(LightCone::const_iterator h=
	lightcone->begin(); h != lightcone->end(); ++h) {
    
    if(h->z < z_min || h->z >= z_max)
      continue;

    hod->compute_param_z(h->z);

    // centrals
    double ncen= hod->ncen(h->M);
    

    if(rand_uniform() > ncen)
      continue;

    p.x[0]= h->x[0];
    p.x[1]= h->x[1];
    p.x[2]= h->x[2];
    p.z   = h->z;
    p.radec[0] = h->radec[0];
    p.radec[1] = h->radec[1];

    p.M   = h->M;
    p.rsat = 0.0f;
    p.vsat = 0.0f;
    p.vr  = h->vr;

    cat->push_back(p);
    ncen_total++;

    double nsat_mean= hod->nsat(h->M);
    int nsat= rand_poisson(nsat_mean);
    //fprintf(stderr, "nsat %f %d\n", nsat_mean, nsat);
    nsat_total += nsat;

    // satellites
    for(int isat=0; isat<nsat; ++isat) {
      satellite(&*h, &p);
      p.rsat = util::norm(p.x);
      p.vsat =  p.vr;
      p.x[0] += h->x[0];
      p.x[1] += h->x[1];
      p.x[2] += h->x[2];
      p.z     = h->z;
      p.vr   += h->vr;
      p.radec[0] = h->radec[0];
      p.radec[1] = h->radec[1];

      cat->push_back(p);
    }    
  }
  cat->ncen= ncen_total;
  cat->nsat= nsat_total;

  msg_printf(msg_verbose,
	     "%d centrals, %d satellites\n", ncen_total, nsat_total);
}


void catalogue_generate_centrals(Hod* const hod,
      LightCone const * const lightcone, const double z_min, const double z_max,
      Catalogue * const cat)
{
  catalogue_init();
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

    if(rand_uniform() > ncen)
      continue;

    p.x[0]= h->x[0];
    p.x[1]= h->x[1];
    p.x[2]= h->x[2];
    p.rsat = 0.0f;
    p.vsat = 0.0f;

    p.M   = h->M;
    p.vr  = h->vr;
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
  catalogue_init();
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
    
    if(rand_uniform() <= ncen) {
      
      p.x[0]= h->x[0];
      p.x[1]= h->x[1];
      p.x[2]= h->x[2];
      p.M= h->M;
      p.rsat = 0.0f;
      p.vsat = 0.0f;
      p.vr  = h->vr;
      p.z   = h->z;
      p.radec[0] = h->radec[0];
      p.radec[1] = h->radec[1];

      cat->push_back(p);
      ncen_total++;

      double nsat_mean= hod->nsat(h->M);
      int nsat= rand_poisson(nsat_mean);

      if(nsat > 0) {
	qsat[iz].push(h->M, h->rs, nsat);
	nsat_total += nsat;
      }
    }
    else if(!qsat[iz].empty()) {
      // pop satellites
      Halo hh= *h;
      qsat[iz].pop(hh);
      
      satellite(&hh, &p);
      p.rsat = util::norm(p.x);
      p.vsat = p.vr;
      p.x[0] += h->x[0];
      p.x[1] += h->x[1];
      p.x[2] += h->x[2];
      p.z     = h->z;
      p.vr   += h->vr;
      p.radec[0] = h->radec[0];
      p.radec[1] = h->radec[1];

      cat->push_back(p);
    }
  }
  cat->ncen= ncen_total;
  cat->nsat= nsat_total;
}

double* catalogue_compute_nz(Catalogue const * const cat,
		      const double z_min, const double z_max, const int nbin)
{
  // Return dN/dz
  // Return value: array of dN/nz needs to be freed by the user
  double* count= (double*) calloc(sizeof(double), nbin);
  const double dz= (z_max - z_min)/nbin;
  
  for(Catalogue::const_iterator p= cat->begin(); p != cat->end(); ++p) {
    int i= floor((p->z - z_min)/(z_max - z_min)*nbin);
    if(0 <= i && i < nbin)
      count[i] += 1.0;
  }

  for(int i=0; i<nbin; ++i)
    count[i] /= dz;

  return count;
}

