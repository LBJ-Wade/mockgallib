#include "cola_lightcone.h"
#include "py_assert.h"
#include "py_cola_lightcones.h"


PyObject* py_cola_lightcones_create(PyObject* self, PyObject* args)
{
  // _cola_lightcones_create(_snapshots, _sky, _remap, _slice, _lightcones)
  PyObject *py_snapshots, *py_sky, *py_remap, *py_slice, *py_lightcones;
  
  if(!PyArg_ParseTuple(args, "OOOOO",
		 &py_snapshots, &py_sky, &py_remap, &py_slice, &py_lightcones))
    return NULL;

  Snapshots const * const snapshots= (Snapshots*)
    PyCapsule_GetPointer(py_snapshots, "_Snapshots");
  py_assert_ptr(snapshots);

  LightCones* const lightcones=
    (LightCones*) PyCapsule_GetPointer(py_lightcones, "_LightCones");
  py_assert_ptr(lightcones);

  Sky const * const sky= (Sky*) PyCapsule_GetPointer(py_sky, "_Sky");
  py_assert_ptr(sky);

  Remap const * const remap= (Remap*) PyCapsule_GetPointer(py_remap, "_Remap");
  py_assert_ptr(remap);

  Slice const * const slice= (Slice*) PyCapsule_GetPointer(py_remap, "_Slice");
  py_assert_ptr(slice);

  //Slice slice(remap->boxsize, sky->width, sky->centre);
  
  cola_lightcones_create(snapshots, sky, remap, slice, lightcones);

  Py_RETURN_NONE;
}

