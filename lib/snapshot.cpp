#include <cstdlib>
#include <cstring>
#include "cosmology.h"
#include "snapshot.h"

using namespace std;

Snapshot::Snapshot(const char filename_fof[],
		   const char filename_halo_mass[],
		   const double a_snp_,
		   const double a_min_, const double a_max_) :
  a_snp(a_snp_), a_min(a_min_), a_max(a_max_)
{
  r_min= cosmology_compute_comoving_distance(a_max_);
  r_max= cosmology_compute_comoving_distance(a_min_);

  halo_mass= new HaloMassFoF(filename_halo_mass);

  const int n= strlen(filename_fof);
  filename= (char*) malloc(n+1);
  strncpy(filename, filename_fof, n+1);
}

Snapshot::~Snapshot()
{
  free(filename);
  delete halo_mass;
}

Snapshots::~Snapshots()
{
  for(deque<Snapshot*>::iterator p= begin(); p != end(); ++p)
    delete *p;
}
