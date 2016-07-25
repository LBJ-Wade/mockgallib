//
// wrapping power.cpp
//
#include "Python.h"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"

#include "power.h"
#include "error.h"
#include "py_power.h"
#include "py_assert.h"

using namespace std;

PyMODINIT_FUNC
py_power_module_init()
{
  import_array();

  return NULL;
}


PyObject* py_power_init(PyObject* self, PyObject* args)
{
  //
  // Convert Python string to char* filename
  //
  PyObject* bytes;
  
  if(!PyArg_ParseTuple(args, "O&", PyUnicode_FSConverter, &bytes)) {
    return NULL;
  }
  
  char* filename;
  Py_ssize_t len;
  PyBytes_AsStringAndSize(bytes, &filename, &len);

  try {
    power_init(filename);
  }
  catch(FileNotFoundError) {
    Py_DECREF(bytes);
    PyErr_SetNone(PyExc_FileNotFoundError);
    return NULL;
  }

  Py_DECREF(bytes);

  Py_RETURN_NONE;
}


PyObject* py_power_free(PyObject* self, PyObject* args)
{
  power_free();

  Py_RETURN_NONE;
}


PyObject* py_power_sigma(PyObject* self, PyObject* args)
{
  double R;
  
  if(!PyArg_ParseTuple(args, "d", &R))
    return NULL;

  double sigma= power_compute_sigma(R);

  return Py_BuildValue("d", sigma);
}


PyObject* py_power_n(PyObject* self, PyObject* args)
{
  return Py_BuildValue("i", power_n());
}

  
PyObject* py_power_k(PyObject* self, PyObject* args)
{
  int rank=1;
  npy_intp dims[]= {power_n()};
  PyObject* k = PyArray_SimpleNewFromData(rank, dims, NPY_DOUBLE, power_k());

  return k;
}


PyObject* py_power_P(PyObject* self, PyObject* args)
{
  int rank=1;
  npy_intp dims[]= {power_n()};
  PyObject* P = PyArray_SimpleNewFromData(rank, dims, NPY_DOUBLE, power_P());

  return P;
}


PyObject* py_power_ki(PyObject* self, PyObject* args)
{
  int i;
  
  if(!PyArg_ParseTuple(args, "i", &i))
     return NULL;

  if(0 <= i && i < power_n())
    return Py_BuildValue("d", power_k()[i]);


  PyErr_SetNone(PyExc_IndexError);
  return NULL;
}


PyObject* py_power_Pi(PyObject* self, PyObject* args)
{
  int i;
  
  if(!PyArg_ParseTuple(args, "i", &i))
     return NULL;

  if(0 <= i && i < power_n())
    return Py_BuildValue("d", power_P()[i]);

  PyErr_SetNone(PyExc_IndexError);
  return NULL;
}
