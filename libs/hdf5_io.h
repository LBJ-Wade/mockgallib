#ifndef HDF5_IO_H
#define HDF5_IO_H 1

#include "lightcone.h"

struct HDF5_IO_Params {
  double omega_m;
};

struct LightconeFileError {

};

void hdf5_read_lightcone(const char filename[], LightCone* const p);
//void hdf5_write(const char filename[], const std::vector<Halo>& v, const int slice, HDF5_IO_Params const * const params);

//HDF5_IO_Params* hdf5_read(const char filename[], std::vector<Halo>* const v);

#endif
