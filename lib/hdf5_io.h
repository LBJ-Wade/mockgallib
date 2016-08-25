#ifndef HDF5_IO_H
#define HDF5_IO_H 1

#include "catalogue.h"
#include "lightcone.h"

struct HDF5_IO_Params {
  double omega_m;
};

void hdf5_read_lightcone(const char filename[], LightCone* const p);
void hdf5_write_lightcone(const char filename[], LightCone const * const v);

void hdf5_read_catalogue(const char filename[], Catalogue* const pv);
#endif
