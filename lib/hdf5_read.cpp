#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <hdf5.h>

#include "msg.h"
#include "halo.h"
#include "hdf5_io.h"


using namespace std;

static int read_data_length(hid_t loc, const char name[]);
static void read_data_table(hid_t loc, const char name[],
			    float * const val, const int nx, const int ny,
			    const hsize_t stride, const bool require=true);
static void read_data_scalar(hid_t loc, const char name[], void * const dat,
			     const hid_t mem_type, const hid_t data_type);
static int read_data_int(hid_t loc, const char name[]);
static float read_data_float(hid_t loc, const char name[]);
static double read_data_double(hid_t loc, const char name[]);


void hdf5_read_lightcone(const char filename[], LightCone* const p)
{
  hid_t file= H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  if(file < 0) {
    msg_printf(msg_fatal, "Error: unable to open %s\n", filename);
    throw LightconeFileError();
  }

  H5Eset_auto(H5E_DEFAULT, NULL, NULL);
  
  // Get size of the data
  const int n= read_data_length(file, "x");
  if(n > 0) {
    // prepare vector to read
    Halo halo;
    p->clear();
    p->reserve(n);
    p->insert(p->begin(), n, halo);
    
    assert(!p->empty());
    
    Halo * const h= &p->front();
    
    
    assert(sizeof(Halo) % sizeof(float) == 0);
    const int stride= sizeof(Halo)/sizeof(float);
    
    // positions
    read_data_table(file, "x", h->x, n, 3, stride);
    
    // radial velocities
    read_data_table(file, "vr", &h->vr, n, 1, stride);
    
    // redshift
    read_data_table(file, "z", &h->z, n, 1, stride);
    
    // ra-dec
    read_data_table(file, "ra-dec", h->radec, n, 2, stride);
    
    
    // M200
    read_data_table(file, "M", &h->M, n, 1, stride, false);
    
    // rs
    read_data_table(file, "rs", &h->rs, n, 1, stride, false);
  }
  
  H5Fclose(file);
}

//
// Utilities
//
int read_data_length(hid_t loc, const char name[])
{
  hid_t dataset= H5Dopen(loc, name, H5P_DEFAULT);
  if(dataset < 0) {
    msg_printf(msg_fatal, "Error: unable to open dataset: %s\n", name);
    throw LightconeFileError();
  }

  hid_t dataspace = H5Dget_space(dataset); assert(dataspace > 0);
  
  hsize_t dims[2];
  H5Sget_simple_extent_dims(dataspace, dims, NULL);
  H5Sclose(dataspace);
  H5Dclose(dataset);  

  return dims[0];
}
  
void read_data_table(hid_t loc, const char name[],
		     float * const val, const int nx, const int ny,
		     const hsize_t stride, const bool require)
{
  hid_t dataset= H5Dopen(loc, name, H5P_DEFAULT);

  if(dataset < 0) {
    if(require) {
      msg_printf(msg_fatal, "Error: unable to open dataset: %s\n", name);
      throw LightconeFileError();
    }
    else {
      msg_printf(msg_verbose, "Unable to open dataset %s, but this is not required\n", name);
      return;
    }
  }

  const hsize_t data_size_mem= stride*nx;
  hid_t dataspace_mem= H5Screate_simple(1, &data_size_mem, 0);
  const hsize_t offset= 0;
  const hsize_t block_size= ny;
  const hsize_t block_count= nx;

  H5Sselect_hyperslab(dataspace_mem, H5S_SELECT_SET,
		      &offset, &stride, &block_count, &block_size);


  const herr_t status_read= 
    H5Dread(dataset, H5T_NATIVE_FLOAT, dataspace_mem, H5S_ALL, H5P_DEFAULT, val);

  if(status_read < 0) {
    msg_printf(msg_fatal, "Error: unable to read dataset: %s\n", name);
    throw LightconeFileError();
  }

  H5Sclose(dataspace_mem);
  H5Dclose(dataset);
}

void read_data_scalar(hid_t loc, const char name[], void * const dat,
		       const hid_t mem_type, const hid_t data_type)
{
  const hid_t scalar= H5Screate(H5S_SCALAR);
  hid_t data= H5Dopen(loc, name, H5P_DEFAULT);
  if(data < 0) {
    msg_printf(msg_fatal, "Error: unable to read data: %s\n", name);
    throw LightconeFileError();
  }

  herr_t status= H5Dread(data, mem_type, scalar, H5S_ALL,
			 H5P_DEFAULT, dat);
  assert(status >= 0);

  H5Dclose(data);
  H5Sclose(scalar);
}

int read_data_int(hid_t loc, const char name[])
{
  int dat;
  read_data_scalar(loc, name, &dat, H5T_NATIVE_INT, H5T_STD_I32LE);

  return dat;
}

float read_data_float(hid_t loc, const char name[])
{
  float dat;
  read_data_scalar(loc, name, &dat, H5T_NATIVE_FLOAT, H5T_IEEE_F32LE);

  return dat;
}

double read_data_double(hid_t loc, const char name[])
{
  double dat;
  read_data_scalar(loc, name, &dat, H5T_NATIVE_DOUBLE, H5T_IEEE_F64LE);

  return dat;
}
