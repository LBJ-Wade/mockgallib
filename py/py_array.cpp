//
// Read a ascii file and broadcast to all MPI nodes
//

#include "Python.h"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"

#include <iostream>

#include "msg.h"
#include "comm.h"
#include "error.h"
#include "py_array.h"

using namespace std;

PyMODINIT_FUNC
py_array_module_init()
{
  import_array();

  return NULL;
}

static double* read_file(const char filename[], int* const nrow_out, int* const ncol_out);

PyObject* py_array_loadtxt(PyObject* self, PyObject* args)
{
  // _array_loadtxt(filename)

  int ncol=0, nrow=0;
  double* buf= 0;
  
  if(comm_this_rank() == 0) {
    PyObject* bytes;
    
    if(!PyArg_ParseTuple(args, "O&", PyUnicode_FSConverter, &bytes)) {
      return NULL;
    }
    
    char* filename;
    Py_ssize_t len;
    PyBytes_AsStringAndSize(bytes, &filename, &len);
    
    try {
      buf= read_file(filename, &ncol, &nrow);
    }
    catch(FileNotFoundError) {
      Py_DECREF(bytes);
      PyErr_SetNone(PyExc_FileNotFoundError);
      return NULL;
    }
    catch(IOError) {
      Py_DECREF(bytes);
      PyErr_SetNone(PyExc_IOError);
      return NULL;
    }

    Py_DECREF(bytes);
  }

  // Broadcast data to all MPI nodes
  ncol= comm_bcast_int(ncol);
  nrow= comm_bcast_int(nrow);
  const int n= ncol*nrow;

  npy_intp dims[]= {nrow, ncol};

  PyObject* py_arr= PyArray_SimpleNew(2, dims, NPY_DOUBLE);
  double* arr= (double*) PyArray_DATA((PyArrayObject*) py_arr);

  if(comm_this_rank() == 0) {
    for(int i=0; i<n; ++i)
      arr[i]= buf[i];
    free(buf);
  }

  comm_mpi_bcast_double(arr, n);

  return py_arr;
}

//
// Read space/tab separated ascii file
//
// * line starting with # is a comment
// * All lines not comment mush have same number of values (ncol)
//
static char* skip_space(char* p)
{
  while(*p == ' ' || *p == '\t')
    p++;
  return p;
}

static char* get_next_space(char* p)
{
  while(!(*p == ' ' || *p == '\t' || *p == '\0'))
    p++;

  return p;
}
    
double* read_file(const char filename[], int* const nrow_out, int* const ncol_out)
{
  FILE* const fp= fopen(filename, "r");
  if(fp == 0) {
    msg_printf(msg_fatal, "Unable to open: %s", filename);
    throw FileNotFoundError();
  }

  char line[256];
  int nalloc= 1000;

  int ibuf= 0;
  double* buf= (double*) malloc(sizeof(double)*nalloc);
  int ncol= 0;
  int nrow= 0;
  int ncol_check= 0;
  
  while(fgets(line, 255, fp)) {
    if(line[0] == '#') continue;
    
    line[255]= '\0';
    char* p= line;

    nrow++;
    ncol_check = 0;
    while(1) {
      p= skip_space(p);
      if(!isdigit(*p))
	break;
      
      double a= atof(p);
      if(ibuf == nalloc) {
	msg_printf(msg_debug, "Reallocating array buffer %d -> %d\n",
		 nalloc, 2*nalloc);
	nalloc *= 2;
	buf= (double*) realloc(buf, sizeof(double)*2*nalloc); assert(buf);
      }

      buf[ibuf++]= a;

      if(nrow == 1)
	ncol++;
      else
	ncol_check++;
      
      p= get_next_space(p);
    }

    if(nrow > 1 && ncol != ncol_check) {
      msg_printf(msg_fatal, "Number of columns are not the same\n");
      fclose(fp);
      throw IOError();
    }
  }
    
  fclose(fp);

  *nrow_out= nrow;
  *ncol_out= ncol;
  return buf;
}
