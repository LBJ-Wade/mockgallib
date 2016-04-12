#include "Python.h"
#include "cosmology.h"
#include "const.h"
#include "mf.h"

static void py_mf_free(PyObject *obj);

PyObject* py_mf_alloc(PyObject* self, PyObject* args)
{
  // py_mf_alloc()
  MF* const mf= mf_alloc();

  return PyCapsule_New(mf, "_MassFunction", py_mf_free);
}

void py_mf_free(PyObject *obj)
{
  MF* const mf= (MF*) PyCapsule_GetPointer(obj, "_MassFunction");

  mf_free(mf);
}

PyObject* py_mf_set_redshift(PyObject* self, PyObject* args)
{
  // _mf_set_redshift(_mf, a)
  PyObject* py_mf;
  double a;
  if(!PyArg_ParseTuple(args, "Od", &py_mf, &a))
    return NULL;

  MF* const mf=
    (MF*) PyCapsule_GetPointer(py_mf, "_MassFunction");

  mf_set_redshift(mf, a);

  Py_RETURN_NONE;
}  
  
PyObject* py_mf_f(PyObject* self, PyObject* args)
{
  PyObject* py_mf;
  double nu;
  if(!PyArg_ParseTuple(args, "Od", &py_mf, &nu))
    return NULL;

  MF* const mf=
    (MF*) PyCapsule_GetPointer(py_mf, "_MassFunction");

  double f= mf_f(mf, nu);
  return Py_BuildValue("d", f);
}

