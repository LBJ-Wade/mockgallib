#ifndef SNAPSHOT_H
#define SNAPSHOT_H 1

#include <deque>
#include "halo_mass.h"
#include "mf_cumulative.h"

class Snapshot {
 public:
  Snapshot(const char filename_fof[],
	   const char filename_part[],
	   const char filename_halo_mass[],
	   const double M_part_min, const double M_halo_min, 
	   const double a_snp, const double a_min, const double a_max);
  ~Snapshot();

  void load_halo_mass(const char filename);

  const double a_snp, a_min, a_max;
  const double M_part_min, M_halo_min;
  float r_min, r_max;
  HaloMassFoF const * halo_mass;
  MfCumulative const * mfc;
  char* filename_fof;
  char* filename_part;
};

class Snapshots: public std::deque<Snapshot*>
{
 public:
  ~Snapshots();
};

#endif
