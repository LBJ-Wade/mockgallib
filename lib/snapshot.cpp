#include <cstdlib>
#include <cstring>
#include "cosmology.h"
#include "snapshot.h"

using namespace std;

Snapshot::Snapshot(const char filename_fof_[],
		   const char filename_part_[],
		   const char filename_halo_mass[],
		   const double M_part_min_, const double M_halo_min_,
		   const double a_snp_,
		   const double a_min_, const double a_max_) :
  M_part_min(M_part_min_), M_halo_min(M_halo_min_),
  a_snp(a_snp_), a_min(a_min_), a_max(a_max_)
{
  r_min= cosmology_compute_comoving_distance(a_max_);
  r_max= cosmology_compute_comoving_distance(a_min_);

  halo_mass= new HaloMassFoF(filename_halo_mass);
  mfc= new MfCumulative(a_snp_);

  int n= strlen(filename_fof_);
  filename_fof= (char*) malloc(n+1);
  strncpy(filename_fof, filename_fof_, n+1);

  n= strlen(filename_part_);
  filename_part= (char*) malloc(n+1);
  strncpy(filename_part, filename_part_, n+1);
}

Snapshot::~Snapshot()
{
  free(filename_fof);
  free(filename_part);
  delete halo_mass;
  delete mfc;
}

Snapshots::~Snapshots()
{
  for(deque<Snapshot*>::iterator p= begin(); p != end(); ++p)
    delete *p;
}
