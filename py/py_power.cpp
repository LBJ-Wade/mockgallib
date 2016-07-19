//
// wrapping power.cpp
//
#include "Python.h"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"

#include "power.h"
#include "py_power.h"
#include "py_assert.h"

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
  
  if(!PyArg_ParseTuple(args, "O&", PyUnicode_FSConverter, &bytes)) {
    return NULL;
  }
  
  char* filename;
  Py_ssize_t len;
  PyBytes_AsStringAndSize(bytes, &filename, &len);

  PowerSpectrum* ps;

  try {
    ps= new PowerSpectrum(filename);
  }
  catch(PowerFileError) {
    Py_DECREF(bytes);
    PyErr_SetNone(PyExc_IOError);
    return NULL;
  }

  Py_DECREF(bytes);

  return PyCapsule_New(ps, "_PowerSpectrum", py_power_free);
}


void py_power_free(PyObject *obj)
{
  PowerSpectrum* const ps=
    (PowerSpectrum*) PyCapsule_GetPointer(obj, "_PowerSpectrum");
  py_assert_void(ps);

  delete ps;
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

  double sigma= ps->compute_sigma(R);

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


  PyErr_SetNone(PyExc_IndexError);
  return NULL;
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

  PyErr_SetNone(PyExc_IndexError);
  return NULL;
}
