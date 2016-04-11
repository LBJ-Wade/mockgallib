//
// wrapping power.cpp
//
//#include <iostream>
#include "Python.h"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"

#include "power.h"
#include "py_power.h"

using namespace std;

static void py_power_free(PyObject *obj);

PyMODINIT_FUNC
py_power_module_init()
{
  import_array();

  return NULL;
}

PyObject* py_power_alloc(PyObject* self, PyObject* args)
{
  //
  // Convert Python string to char* filename
  //
  PyObject* bytes;
  char* filename;
  Py_ssize_t len;

  if(!PyArg_ParseTuple(args, "O&", PyUnicode_FSConverter, &bytes)) {
    return NULL;
  }

  PyBytes_AsStringAndSize(bytes, &filename, &len);

  printf("filename= %s\n", filename);

  PowerSpectrum* ps;

  try {
    ps= power_alloc(filename);
  }
  catch(PowerFileError) {
    Py_DECREF(bytes);
    Py_RETURN_NONE;
  }

  Py_DECREF(bytes);

  return PyCapsule_New(ps, "_PowerSpectrum", py_power_free);
}


void py_power_free(PyObject *obj)
{
  PowerSpectrum* const ps=
    (PowerSpectrum*) PyCapsule_GetPointer(obj, "_PowerSpectrum");

  power_free(ps);
}


PyObject* py_power_sigma(PyObject* self, PyObject* args)
{
  PyObject* py_ps;
  double R;
  
  if(!PyArg_ParseTuple(args, "Od", &py_ps, &R))
    return NULL;

  PowerSpectrum* ps;
  if (!(ps =  (PowerSpectrum *) PyCapsule_GetPointer(py_ps, "_PowerSpectrum")))
    return NULL;

  double sigma= power_sigma(ps, R);

  return Py_BuildValue("d", sigma);
}

PyObject* py_power_n(PyObject* self, PyObject* args)
{
  PyObject* py_ps;
  
  if(!PyArg_ParseTuple(args, "O", &py_ps))
     return NULL;

  PowerSpectrum* ps;
  if (!(ps =  (PowerSpectrum *) PyCapsule_GetPointer(py_ps, "_PowerSpectrum")))
    return NULL;

  return Py_BuildValue("i", ps->n);
}

  
PyObject* py_power_k(PyObject* self, PyObject* args)
{
  PyObject* py_ps;
  
  if(!PyArg_ParseTuple(args, "O", &py_ps))
     return NULL;

  PowerSpectrum* ps;
  if (!(ps =  (PowerSpectrum *) PyCapsule_GetPointer(py_ps, "_PowerSpectrum")))
    return NULL;

  int rank=1;
  npy_intp dims[]= {ps->n};
  
  PyObject* k = PyArray_SimpleNewFromData(rank, dims, NPY_DOUBLE, ps->k);
  return k;
}

PyObject* py_power_P(PyObject* self, PyObject* args)
{
  PyObject* py_ps;
  
  if(!PyArg_ParseTuple(args, "O", &py_ps))
     return NULL;

  PowerSpectrum* ps;
  if (!(ps =  (PowerSpectrum *) PyCapsule_GetPointer(py_ps, "_PowerSpectrum")))
    return NULL;

  int rank=1;
  npy_intp dims[]= {ps->n};
  
  PyObject* P = PyArray_SimpleNewFromData(rank, dims, NPY_DOUBLE, ps->P);
  return P;
}

PyObject* py_power_ki(PyObject* self, PyObject* args)
{
  PyObject* py_ps;
  int i;
  
  if(!PyArg_ParseTuple(args, "Oi", &py_ps, &i))
     return NULL;

  PowerSpectrum* ps;
  if (!(ps =  (PowerSpectrum *) PyCapsule_GetPointer(py_ps, "_PowerSpectrum")))
    return NULL;

  if(0 <= i && i < ps->n)
    return Py_BuildValue("d", ps->k[i]);

  Py_RETURN_NONE;
}

PyObject* py_power_Pi(PyObject* self, PyObject* args)
{
  PyObject* py_ps;
  int i;
  
  if(!PyArg_ParseTuple(args, "Oi", &py_ps, &i))
     return NULL;

  PowerSpectrum* ps;
  if (!(ps =  (PowerSpectrum *) PyCapsule_GetPointer(py_ps, "_PowerSpectrum")))
    return NULL;

  if(0 <= i && i < ps->n)
    return Py_BuildValue("d", ps->P[i]);

  Py_RETURN_NONE;
}

/*
static PyMethodDef Methods[] = {
  {"_power_alloc", py_power_alloc, METH_VARARGS,
   "read and allocate power spectrum"},
  {"_power_sigma", py_power_sigma, METH_VARARGS,
   "compute sigma(R); density fluctuation smoothed with radius R"},
  {"_power_n", py_power_n, METH_VARARGS,
   "return number of power spectrum rows"},
  {"_power_k", py_power_k, METH_VARARGS,
   "return wavenumber array k"},
  {"_power_P", py_power_P, METH_VARARGS,
   "return power spectrum array P"},
  {"_power_ki", py_power_ki, METH_VARARGS,
   "return k[i]"},
  {"_power_Pi", py_power_Pi, METH_VARARGS,
   "return P[i]"},
  {NULL, NULL, 0, NULL}
};

void py_power_add_methods(vector<PyMethodDef>& methods)
{
  int i=0;
  while(Methods[i].ml_name != NULL) {
    methods.push_back(Methods[i++]);
  }
}


static struct PyModuleDef module = {
  PyModuleDef_HEAD_INIT,
  "mockgallib.powerspectrum", // name of this module
  "power spectrum module", // Doc String
  -1,
  Methods
};


PyMODINIT_FUNC
PyInit__mockgallib_power(void) {
  cerr << "py_power initialised\n";
  import_array();
  return PyModule_Create(&module);
}

*/
