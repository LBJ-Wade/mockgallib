#include "nbar.h"
#include "py_nbar.h"

static void py_nbar_free(PyObject *obj);

PyObject* py_nbar_alloc(PyObject* self, PyObject* args)
{
  // _nbar_alloc(_ps, _hod)
  
  PyObject *py_ps, *py_hod;
  if(!PyArg_ParseTuple(args, "OO", &py_ps, &py_hod))
    return NULL;

  PowerSpectrum* const ps=
    (PowerSpectrum*) PyCapsule_GetPointer(py_ps, "_PowerSpectrum");

  Hod* const hod=
    (Hod*) PyCapsule_GetPointer(py_hod, "_HOD");

  NbarIntegration* const ni= nbar_integration_alloc(ps, hod);

  //return PyCapsule_New(ni, "_NbarIntegration", ps);
  return PyCapsule_New(ni, "_NbarIntegration", py_nbar_free);
  //Py_RETURN_NONE; // debug!!!
}

void py_nbar_free(PyObject *obj)
{
  NbarIntegration* const ni=
      (NbarIntegration*) PyCapsule_GetPointer(obj, "_NbarIntegration");

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
