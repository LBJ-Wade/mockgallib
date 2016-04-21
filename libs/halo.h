#ifndef HALO_H
#define HALO_H 1

struct Halo {
  int nfof;
  float x[3], v[3], r;
  float M;
  float ra, dec;
  int slice;
  float a; // scale factor of the snapshot
  float rs; // NFW rs, physical 1/h kpc
  float z; // redshift at radius r
};

#endif
