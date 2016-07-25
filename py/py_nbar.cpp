#include "msg.h"
#include "nbar.h"
#include "py_nbar.h"

static void py_nbar_free(PyObject *obj);

PyObject* py_nbar_alloc(PyObject* self, PyObject* args)
{
  // _nbar_alloc(_ps, _hod)
  
  PyObject *py_hod;
  if(!PyArg_ParseTuple(args, "O", &py_hod))
    return NULL;

  Hod* const hod=
    (Hod*) PyCapsule_GetPointer(py_hod, "_HOD");

  NbarIntegration* const ni= nbar_integration_alloc(hod);

  return PyCapsule_New(ni, "_NbarIntegration", py_nbar_free);
}


void py_nbar_free(PyObject *obj)
{
  NbarIntegration* const ni=
      (NbarIntegration*) PyCapsule_GetPointer(obj, "_NbarIntegration");

  msg_printf(msg_debug, "Freeing nbar\n");
  
  nbar_integration_free(ni);
}

PyObject* py_nbar_compute(PyObject* self, PyObject* args)
{
  // _nbar_compute(_ni)

  PyObject* py_ni;
  double z;
  if(!PyArg_ParseTuple(args, "Od", &py_ni, &z))
    return NULL;

  NbarIntegration* const ni=
    (NbarIntegration*) PyCapsule_GetPointer(py_ni, "_NbarIntegration");

  double nbar= nbar_compute(ni, z);

  return Py_BuildValue("d", nbar);
}
