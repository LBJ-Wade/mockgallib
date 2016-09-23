#include "mf_cumulative.h"
#include "py_mf_cumulative.h"
#include "py_assert.h"

static void py_mf_cumulative_free(PyObject *obj);

PyObject* py_mf_cumulative_alloc(PyObject* self, PyObject* args)
{
  // py_mf_cumulatice_alloc(a)
  double a;
  if(!PyArg_ParseTuple(args, "d", &a)) {
    return NULL;
  }

  MfCumulative* const mfc= new MfCumulative(a);

  return PyCapsule_New(mfc, "_MfCumulatice", py_mf_cumulative_free);
}

void py_mf_cumulative_free(PyObject *obj)
{
  MfCumulative* const mfc=
    (MfCumulative*) PyCapsule_GetPointer(obj, "_MfCumulative");
  py_assert_void(mfc);

  delete mfc;
}

