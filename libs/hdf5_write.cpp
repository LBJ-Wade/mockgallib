#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <hdf5.h>

#include "halo.h"
#include "hdf5_io.h"


using namespace std;

static void write_data_int(hid_t loc, const char name[], const int val);
static void write_data_float(hid_t loc, const char name[], const float val);
static void write_data_double(hid_t loc, const char name[], const double val);
static void write_data_table(hid_t loc, const char name[],
		     float const * const val, const int nx, const int ny);

static inline float norm(const float x[])
{
  return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}

static inline float dot(const float x[], const float y[])
{
  return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

//typedef float[3] float3;


void hdf5_write(const char filename[], const vector<Halo>& v, const int slice, HDF5_IO_Params const * const params)
{
  // slice: slice number to write. Writes all slices if slice<0
  
  hid_t file= H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if(file < 0) {
    cerr << "Error: unable to create: " << filename << endl;
    throw filename;
  }

  // parameters
  hid_t group= H5Gcreate(file, "parameters",
			 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // comment: good to have z range, too?
  write_data_double(group, "omega_m", params->omega_m);
  
  
  vector<float> data;

  // positions
  int n= 0;
  for(vector<Halo>::const_iterator h= v.begin(); h != v.end(); ++h) {
    if(slice >= 0 && h->slice != slice) continue;

    data.push_back(h->x[0]);
    data.push_back(h->x[1]);
    data.push_back(h->x[2]);
    n++;
  }
  assert(data.size() == 3*n);
  write_data_table(file, "x", &data.front(), n, 3); data.clear();

  // radial velocities
  for(vector<Halo>::const_iterator h= v.begin(); h != v.end(); ++h) {
    if(slice >= 0 && h->slice != slice) continue;

    data.push_back(dot(h->x, h->v)/norm(h->x));
  }
  assert(data.size() == n);
  write_data_table(file, "vr", &data.front(), n, 1); data.clear();

  // redshift
  for(vector<Halo>::const_iterator h= v.begin(); h != v.end(); ++h) {
    if(slice >= 0 && h->slice != slice) continue;

    data.push_back(h->z);
  }
  assert(data.size() == n);
  write_data_table(file, "z", &data.front(), n, 1); data.clear();

  // ra-dec
  for(vector<Halo>::const_iterator h= v.begin(); h != v.end(); ++h) {
    if(slice >= 0 && h->slice != slice) continue;

    data.push_back(h->ra);
    data.push_back(h->dec);
  }
  assert(data.size() == 2*n);  
  write_data_table(file, "ra-dec", &data.front(), n, 2); data.clear();


  if(v.front().M > 0.0f) {
    // M200
    for(vector<Halo>::const_iterator h= v.begin(); h != v.end(); ++h) {
      if(slice >= 0 && h->slice != slice) continue;
      
      data.push_back(h->M);
    }
    assert(data.size() == n);
    write_data_table(file, "M", &data.front(), n, 1); data.clear();

    // rs
    for(vector<Halo>::const_iterator h= v.begin(); h != v.end(); ++h) {
      if(slice >= 0 && h->slice != slice) continue;
      
      data.push_back(h->rs);
    }
    assert(data.size() == n);
    write_data_table(file, "rs", &data.front(), n, 1); data.clear();
  }
  
  //write_data_float(group, "omega_m", omega_m);
  //write_data_float(group, "lambda", lambda);


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
    cerr << "Error: unable to create int data: " <<  name << endl;
    throw name;
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
    cerr << "Error: unable to create float data: " << name << endl;
    throw name;
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
    cerr << "Error: unable to create float data: " << name << endl;
    throw name;
  }

  herr_t status= H5Dwrite(data, H5T_NATIVE_DOUBLE, scalar, H5S_ALL,
			  H5P_DEFAULT, &val);
  assert(status >= 0);

  H5Dclose(data);
  H5Sclose(scalar);
}

void write_data_table(hid_t loc, const char name[],
		      float const * const val, const int nx, const int ny)
{
  const hsize_t rank= 2;
  const hsize_t data_size_file[]= {nx, ny};

  hid_t dataspace_file= H5Screate_simple(rank, data_size_file, 0);
  hid_t dataset= H5Dcreate(loc, name, H5T_IEEE_F32LE, dataspace_file,
			   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  if(dataset < 0) {
    cerr << "Error: unable to create dataset: " << name << endl;
    throw name;
  }

  const hsize_t stride= ny;
  const hsize_t data_size_mem= nx*ny;
  hid_t dataspace_mem= H5Screate_simple(1, &data_size_mem, 0);
  const hsize_t offset= 0;
  const hsize_t block_size= ny;
  const hsize_t block_count= nx;

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
