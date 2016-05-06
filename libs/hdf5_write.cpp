#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <hdf5.h>

#include "msg.h"
#include "halo.h"
#include "hdf5_io.h"


using namespace std;

static void write_data_int(hid_t loc, const char name[], const int val);
static void write_data_float(hid_t loc, const char name[], const float val);
static void write_data_double(hid_t loc, const char name[], const double val);
static void write_data_table(hid_t loc, const char name[],
		float const * const val,
		const int nrow, const int ncol, const hsize_t stride);

static inline float norm(const float x[])
{
  return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}

static inline float dot(const float x[], const float y[])
{
  return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

//typedef float[3] float3;

void hdf5_write_lightcone(const char filename[], LightCone const * const v)
{
  hid_t file= H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if(file < 0) {
    msg_printf(msg_error, "Error: unable to create: %s\n", filename);
    throw LightconeFileError();
  }

  Halo const * const h= &v->front();

  assert(sizeof(Halo) % sizeof(float) == 0);
  const int stride= sizeof(Halo) / sizeof(float);
  const int n= v->size();
  
  // positions
  write_data_table(file, "x", h->x, n, 3, stride);

  // radial velocities
  write_data_table(file, "vr", &h->vr, n, 1, stride);

  // redshift
  write_data_table(file, "z", &h->z, n, 1, stride);

  // ra-dec
  write_data_table(file, "ra-dec", h->radec, n, 2, stride);

  if(h->M > 0.0f) {
    // M200
    write_data_table(file, "M", &h->M, n, 1, stride);

    // rs
    write_data_table(file, "rs", &h->rs, n, 1, stride);
  }
  
  H5Fclose(file);
}

//
// Utilities
//
void write_data_int(hid_t loc, const char name[], const int val)
{
  const hid_t scalar= H5Screate(H5S_SCALAR);
  hid_t data= H5Dcreate(loc, name, H5T_STD_I32LE, scalar, 
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if(data < 0) {
    msg_printf(msg_error, "Error: unable to create int data: %s\n", name);
    throw LightconeFileError();
  }

  herr_t status= H5Dwrite(data, H5T_NATIVE_INT, scalar, H5S_ALL,
			  H5P_DEFAULT, &val);
  assert(status >= 0);

  H5Dclose(data);
  H5Sclose(scalar);
}

void write_data_float(hid_t loc, const char name[], const float val)
{
  const hid_t scalar= H5Screate(H5S_SCALAR);
  hid_t data= H5Dcreate(loc, name, H5T_IEEE_F32LE, scalar, 
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if(data < 0) {
    msg_printf(msg_error, "Error: unable to create float data: %s\n", name);
    throw LightconeFileError();
  }

  herr_t status= H5Dwrite(data, H5T_NATIVE_FLOAT, scalar, H5S_ALL,
			  H5P_DEFAULT, &val);
  assert(status >= 0);

  H5Dclose(data);
  H5Sclose(scalar);
}

void write_data_double(hid_t loc, const char name[], const double val)
{
  const hid_t scalar= H5Screate(H5S_SCALAR);
  hid_t data= H5Dcreate(loc, name, H5T_IEEE_F64LE, scalar, 
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if(data < 0) {
    msg_printf(msg_error, "Error: unable to create float data: %s\n", name);
    throw LightconeFileError();
  }

  herr_t status= H5Dwrite(data, H5T_NATIVE_DOUBLE, scalar, H5S_ALL,
			  H5P_DEFAULT, &val);
  assert(status >= 0);

  H5Dclose(data);
  H5Sclose(scalar);
}

void write_data_table(hid_t loc, const char name[], float const * const val,
		      const int nrow, const int ncol, const hsize_t stride)
{
  const hsize_t rank= ncol == 1 ? 1 : 2;
  const hsize_t data_size_file[]= {nrow, ncol};

  hid_t dataspace_file= H5Screate_simple(rank, data_size_file, 0);
  hid_t dataset= H5Dcreate(loc, name, H5T_IEEE_F32LE, dataspace_file,
			   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  if(dataset < 0) {
    msg_printf(msg_error, "Error: unable to create dataset: %s\n", name);
    throw LightconeFileError();
  }

  const hsize_t data_size_mem= nrow*stride;
  hid_t dataspace_mem= H5Screate_simple(1, &data_size_mem, 0);
  const hsize_t offset= 0;
  const hsize_t block_size= ncol;
  const hsize_t block_count= nrow;

  H5Sselect_hyperslab(dataspace_mem, H5S_SELECT_SET,
		      &offset, &stride, &block_count, &block_size);

    
  const herr_t status_write= H5Dwrite(dataset, H5T_NATIVE_FLOAT, 
				      dataspace_mem, dataspace_file,
				      H5P_DEFAULT, val);

  assert(status_write >= 0);
  H5Sclose(dataspace_mem);
  H5Dclose(dataset);
  H5Sclose(dataspace_file);
}
