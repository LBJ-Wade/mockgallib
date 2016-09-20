#ifndef MF_CUMULATIVE
#define MF_CUMULATIVE 1

#include "Python.h"
#include "mf_cumulative.h"
#include "py_assert.h"

static void py_mf_cumulative_free(PyObject *obj);

PyObject* py_mf_cumulative_alloc(PyObject* self, PyObject* args)
{
  // py_mf_alloc(a)
  double a;
  if(!PyArg_ParseTuple(args, "d", &a)) {
    return NULL;
  }

  MfCumulative* const mfc= new MfCumulative(a);

  return PyCapsule_New(mfc, "_MfCumulative", py_mf_cumulative_free);
}

void py_mf_cumulative_free(PyObject *obj)
{
  MF* const mfc= (MF*) PyCapsule_GetPointer(obj, "_MfCumulative");
  py_assert_void(mfc);

  delete mfc;
}

PyObject* py_mf_cumulative_M(PyObject* self, PyObject* args)
{
  // py_mf_cumulative_M(_mfc, n)
  // returns M: n(>M) = n
  PyObject* py_mfc;
  double n;
  if(!PyArg_ParseTuple(args, "Od", &py_mfc, &n)) {
    return NULL;
  }

  MfCumulative* const mfc=
    (MfCumulative*) PyCapsule_GetPointer(py_mfc, "_MfCumulative");
  py_assert_ptr(mfc);


  return Py_BuildValue("d", mfc->M(n));
}




#endif
