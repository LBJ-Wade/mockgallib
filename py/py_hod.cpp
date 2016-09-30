#include "hod.h"
#include "py_hod.h"
#include "py_assert.h"

static void py_hod_free(PyObject *obj);

PyObject* py_hod_alloc(PyObject* self, PyObject* args)
{
  double z0;
  if(!PyArg_ParseTuple(args, "d", &z0))
    return NULL;

  Hod* hod= new Hod(z0); py_assert_ptr(hod)

  return PyCapsule_New(hod, "_HOD", py_hod_free);
}

void py_hod_free(PyObject *obj)
{
  Hod* const hod= (Hod*) PyCapsule_GetPointer(obj, "_HOD");
  py_assert_void(hod);

  delete hod;
}

PyObject* py_hod_get_coef_all(PyObject* self, PyObject* args)
{
  PyObject* py_hod;
  if(!PyArg_ParseTuple(args, "O", &py_hod))
    return NULL;

  Hod* const hod= (Hod*) PyCapsule_GetPointer(py_hod, "_HOD");
  py_assert_ptr(hod);
  
  const int n= hod->n;
  double const * const c= hod->c;

  
  PyObject* const list= PyList_New(n);
  for(int i=0; i<n; ++i)
    PyList_SetItem(list, i, Py_BuildValue("d", c[i]));

  return list;
}

PyObject* py_hod_get_coef(PyObject* self, PyObject* args)
{
  // c[i] =_hod_get_coef(_hod, i)
  PyObject* py_hod;
  int i;

  if(!PyArg_ParseTuple(args, "Oi", &py_hod, &i))
    return NULL;

  Hod* const hod= (Hod*) PyCapsule_GetPointer(py_hod, "_HOD");
  py_assert_ptr(hod);
  
  const int n= hod->n;
  double const * const c= hod->c;

  if(0 <= i && i < n) {
    return Py_BuildValue("d", c[i]);
  }

  PyErr_SetNone(PyExc_IndexError);
  return NULL;
}


PyObject* py_hod_set_coef(PyObject* self, PyObject* args)
{
  // _hod_set_coef(_hod, i, c[i])
  PyObject *py_hod;
  int i;
  double val;
  
  if(!PyArg_ParseTuple(args, "Oid", &py_hod, &i, &val))
    return NULL;

  Hod* const hod= (Hod*) PyCapsule_GetPointer(py_hod, "_HOD");
  py_assert_ptr(hod);
  
  const int n= hod->n;
  double * const c= hod->c;

  if(0 <= i && i < n) {
    c[i]= val;
  }
  else {
    PyErr_SetString(PyExc_ValueError, "HOD coefficient out of range");
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

  Hod* const hod= (Hod*) PyCapsule_GetPointer(py_hod, "_HOD");
  py_assert_ptr(hod);
  
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

  Hod* const hod= (Hod*) PyCapsule_GetPointer(py_hod, "_HOD");
  py_assert_ptr(hod);
  
  hod->compute_param_z(z);
  double nsat= hod->nsat(M);

  return Py_BuildValue("d", nsat);
}

PyObject* py_hod_get_z0(PyObject* self, PyObject* args)
{
  PyObject* py_hod;
  if(!PyArg_ParseTuple(args, "O", &py_hod))
    return NULL;

  Hod* const hod= (Hod*) PyCapsule_GetPointer(py_hod, "_HOD");
  py_assert_ptr(hod);
  
  return Py_BuildValue("d", hod->z0);
}

PyObject* py_hod_compute_logMmin(PyObject* self, PyObject* args)
{
  PyObject* py_hod;
  double z;
  if(!PyArg_ParseTuple(args, "Od", &py_hod, &z))
    return NULL;

  Hod* const hod= (Hod*) PyCapsule_GetPointer(py_hod, "_HOD");
  py_assert_ptr(hod);
  
  return Py_BuildValue("d", hod->compute_logMmin(z));
}

PyObject* py_hod_compute_M1(PyObject* self, PyObject* args)
{
  PyObject* py_hod;
  double z;
  if(!PyArg_ParseTuple(args, "Od", &py_hod, &z))
    return NULL;

  Hod* const hod= (Hod*) PyCapsule_GetPointer(py_hod, "_HOD");
  py_assert_ptr(hod);
  
  return Py_BuildValue("d", hod->compute_M1(z));
}

PyObject* py_hod_compute_sigma(PyObject* self, PyObject* args)
{
  PyObject* py_hod;
  double z;
  if(!PyArg_ParseTuple(args, "Od", &py_hod, &z))
    return NULL;

  Hod* const hod= (Hod*) PyCapsule_GetPointer(py_hod, "_HOD");
  py_assert_ptr(hod);
  
  return Py_BuildValue("d", hod->compute_sigma(z));
}  

PyObject* py_hod_compute_alpha(PyObject* self, PyObject* args)
{
  PyObject* py_hod;
  double z;
  if(!PyArg_ParseTuple(args, "Od", &py_hod, &z))
    return NULL;

  Hod* const hod= (Hod*) PyCapsule_GetPointer(py_hod, "_HOD");
  py_assert_ptr(hod);
  
  return Py_BuildValue("d", hod->compute_alpha(z));
}  
