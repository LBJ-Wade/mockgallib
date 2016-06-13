#ifndef HALO_H
#define HALO_H 1

struct Halo {
  int nfof;
  float x[3], v[3], r, vr;
  float M;
  float radec[2];
  int slice;
  float a; // scale factor of the snapshot
  float rs; // NFW rs, physical 1/h kpc
  float z; // redshift at radius r
};

#endif
