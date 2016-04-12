//
// wrapping sigma.cpp
//
#include <iostream>
#include <cassert>

#include "Python.h"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"

#include "power.h"
#include "sigma.h"
#include "py_sigma.h"

using namespace std;
static void py_sigma_free(PyObject *obj);

PyMODINIT_FUNC
py_sigma_module_init()
{
  import_array();

  return NULL;
}

PyObject* py_sigma_alloc(PyObject* self, PyObject* args)
{
  //
  //
  PyObject* py_ps;
  if(!PyArg_ParseTuple(args, "O", &py_ps))
    return NULL;

  PowerSpectrum* const ps=
    (PowerSpectrum*) PyCapsule_GetPointer(py_ps, "_PowerSpectrum");
  assert(ps);

  Sigma* const s= sigma_alloc(ps);

  return PyCapsule_New(s, "_Sigma", py_sigma_free);
}


void py_sigma_free(PyObject *obj)
{
  Sigma* const s=
    (Sigma*) PyCapsule_GetPointer(obj, "_Sigma");

  sigma_free(s);
}

PyObject* py_sigma_n(PyObject* self, PyObject* args)
{
  PyObject* py_sigma;
  
  if(!PyArg_ParseTuple(args, "O", &py_sigma))
     return NULL;

  Sigma* s;
  if(!(s = (Sigma *) PyCapsule_GetPointer(py_sigma, "_Sigma")))
    return NULL;

  return Py_BuildValue("i", s->n);
}



PyObject* py_sigma_M(PyObject* self, PyObject* args)
{
  // py_sigma_M(_Sigma py_sigma, sigma0) = M(sigma0)
  
  PyObject* py_sigma;
  double sigma0;
  
  if(!PyArg_ParseTuple(args, "Od", &py_sigma, &sigma0))
    return NULL;

  Sigma* const s=
    (Sigma*) PyCapsule_GetPointer(py_sigma, "_Sigma");

  double M= sigma_M(s, sigma0);

  return Py_BuildValue("d", M);
}


PyObject* py_sigma_0inv(PyObject* self, PyObject* args)
{
  // py_sigma_0inv(_Sigma py_sigma, M) = 1.0/sigma0(M)
  
  PyObject* py_sigma;
  double M;
  
  if(!PyArg_ParseTuple(args, "Od", &py_sigma, &M))
    return NULL;

  Sigma* const s=
    (Sigma*) PyCapsule_GetPointer(py_sigma, "_Sigma");

  double sinv= sigma_0inv(s, M);

  return Py_BuildValue("d", sinv);
}

PyObject* py_sigma_M_range(PyObject* self, PyObject* args)
{
  // py_sigma_0inv(_Sigma py_sigma, M) = 1.0/sigma0(M)
  
  PyObject* py_sigma;
  
  if(!PyArg_ParseTuple(args, "O", &py_sigma))
    return NULL;

  Sigma* const s=
    (Sigma*) PyCapsule_GetPointer(py_sigma, "_Sigma");

  return Py_BuildValue("(dd)", s->M_min, s->M_max);
}

PyObject* py_sigma_M_array(PyObject* self, PyObject* args)
{
  PyObject* py_sigma;
  
  if(!PyArg_ParseTuple(args, "O", &py_sigma))
    return NULL;

  Sigma* const s=
    (Sigma*) PyCapsule_GetPointer(py_sigma, "_Sigma");

  int rank=1;
  npy_intp dims[]= {s->n};
  
  PyObject* M = PyArray_SimpleNewFromData(rank, dims, NPY_DOUBLE, s->M);
  return M;
}

PyObject* py_sigma_sinv_array(PyObject* self, PyObject* args)
{
  PyObject* py_sigma;
  
  if(!PyArg_ParseTuple(args, "O", &py_sigma))
    return NULL;

  Sigma* const s=
    (Sigma*) PyCapsule_GetPointer(py_sigma, "_Sigma");

  int rank=1;
  npy_intp dims[]= {s->n};
  
  PyObject* sinv = PyArray_SimpleNewFromData(rank, dims, NPY_DOUBLE, s->sinv);
  return sinv;
}

