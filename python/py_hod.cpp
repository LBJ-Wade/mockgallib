#include <iostream>
#include "hod.h"
#include "py_hod.h"

using namespace std;

static void py_hod_free(PyObject *obj);

PyObject* py_hod_alloc(PyObject* self, PyObject* args)
{
  //
  // Convert Python string to char* filename
  //
  Hod* hod= new Hod();

  return PyCapsule_New(hod, "_HOD", py_hod_free);
}

void py_hod_free(PyObject *obj)
{
  Hod* const hod=
    (Hod*) PyCapsule_GetPointer(obj, "_HOD");

  delete hod;
}

PyObject* py_hod_get_coef(PyObject* self, PyObject* args)
{
  PyObject* py_hod;
  if(!PyArg_ParseTuple(args, "O", &py_hod))
    return NULL;

  Hod* const hod=
    (Hod*) PyCapsule_GetPointer(py_hod, "_HOD");

  const int n= hod->n;
  double const * const c= hod->c;

  
  PyObject* const list= PyList_New(n);
  for(int i=0; i<n; ++i)
    PyList_SetItem(list, i, Py_BuildValue("d", c[i]));

  return list;
}

PyObject* py_hod_set_coef(PyObject* self, PyObject* args)
{
  // _hod_set_coef(_hod, i, c[i])
  PyObject *py_hod;
  int i;
  double val;
  
  if(!PyArg_ParseTuple(args, "Oid", &py_hod, &i, &val))
    return NULL;

  Hod* const hod=
    (Hod*) PyCapsule_GetPointer(py_hod, "_HOD");

  const int n= hod->n;
  double * const c= hod->c;

  cerr << "setting " << i << " " << val << endl;
  
  if(0 <= i && i < n) {
    c[i]= val;
  }

  Py_RETURN_NONE;
}

PyObject* py_hod_ncen(PyObject* self, PyObject* args)
{
  // py_hod_ncen(_hod, M, z)
  
  double M, z;
  PyObject* py_hod;
  if(!PyArg_ParseTuple(args, "Odd", &py_hod, &M, &z))
    return NULL;

  Hod* const hod=
    (Hod*) PyCapsule_GetPointer(py_hod, "_HOD");

  hod->compute_param_z(z);
  double ncen= hod->ncen(M);

  return Py_BuildValue("d", ncen);
}

PyObject* py_hod_nsat(PyObject* self, PyObject* args)
{
  // py_hod_nsat(_hod, M, z)
  
  double M, z;
  PyObject* py_hod;
  if(!PyArg_ParseTuple(args, "Odd", &py_hod, &M, &z))
    return NULL;

  Hod* const hod=
    (Hod*) PyCapsule_GetPointer(py_hod, "_HOD");

  hod->compute_param_z(z);
  double nsat= hod->nsat(M);

  return Py_BuildValue("d", nsat);
}
