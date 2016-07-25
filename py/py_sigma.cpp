//
// wrapping sigma.cpp
//
#include <iostream>
#include <cassert>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"

#include "power.h"
#include "sigma.h"
#include "py_sigma.h"
#include "py_assert.h"

using namespace std;

PyMODINIT_FUNC
py_sigma_module_init()
{
  import_array();

  return NULL;
}


PyObject* py_sigma_init(PyObject* self, PyObject* args)
{
  // _sigma_alloc(_ps, M_min, M_max, n)

  double M_min, M_max;
  int n;
  
  if(!PyArg_ParseTuple(args, "ddi", &M_min, &M_max, &n))
    return NULL;

  Py_RETURN_NONE;
}


PyObject* py_sigma_free(PyObject* self, PyObject* args)
{
  sigma_free();

  Py_RETURN_NONE;
}

PyObject* py_sigma_n(PyObject* self, PyObject* args)
{
  return Py_BuildValue("i", sigma_n());
}


PyObject* py_sigma_M(PyObject* self, PyObject* args)
{
  // py_sigma_M(sigma0) = M(sigma0)
  
  double sigma0;
  
  if(!PyArg_ParseTuple(args, "d", &sigma0))
    return NULL;

  return Py_BuildValue("d", sigma_M(sigma0));
}


PyObject* py_sigma_0inv(PyObject* self, PyObject* args)
{
  // py_sigma_0inv(M) = 1.0/sigma0(M)
  
  double M;
  
  if(!PyArg_ParseTuple(args, "d", &M))
    return NULL;

  return Py_BuildValue("d", sigma_inv(M));
}

PyObject* py_sigma_M_range(PyObject* self, PyObject* args)
{
  // py_sigma_0inv(_Sigma py_sigma, M) = 1.0/sigma0(M)
  
  return Py_BuildValue("(dd)", sigma_M_min(), sigma_M_max());
}

PyObject* py_sigma_M_array(PyObject* self, PyObject* args)
{
  int rank=1;
  npy_intp dims[]= {sigma_n()};
  
  return PyArray_SimpleNewFromData(rank, dims, NPY_DOUBLE, sigma_M_array());
}

PyObject* py_sigma_sinv_array(PyObject* self, PyObject* args)
{
  int rank=1;
  npy_intp dims[]= {sigma_n()};
  
  return PyArray_SimpleNewFromData(rank, dims, NPY_DOUBLE, sigma_sinv_array());
}
