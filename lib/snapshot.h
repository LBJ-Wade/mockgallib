#ifndef SNAPSHOT_H
#define SNAPSHOT_H 1

#include <deque>
#include "halo_mass.h"

class Snapshot {
 public:
  Snapshot(const char filename_fof[],
	   const char filename_halo_mass[],
	   const double a_snp, const double a_min, const double a_max);
  ~Snapshot();

  void load_halo_mass(const char filename);

  const double a_snp, a_min, a_max;
  float r_min, r_max;
  HaloMassFoF const * halo_mass;
  char* filename;
};

class Snapshots: public std::deque<Snapshot*>
{
 public:
  ~Snapshots();
};

#endif
