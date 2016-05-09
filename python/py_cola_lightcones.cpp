#include "cola_file.h"
#include "cola_lightcone.h"
#include "py_assert.h"
#include "py_cola_lightcones.h"


PyObject* py_cola_lightcones_create(PyObject* self, PyObject* args)
{
  // _cola_lightcones_create(_snapshots, _sky, _remap, _slice, _lightcones,
  //                         _random)
  PyObject *py_snapshots, *py_sky, *py_remap, *py_slice, *py_lightcones, *py_random;
  
  if(!PyArg_ParseTuple(args, "OOOOOO",
		       &py_snapshots, &py_sky, &py_remap, &py_slice, &py_lightcones, &py_random))
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

  Slice const * const slice= (Slice*) PyCapsule_GetPointer(py_slice, "_Slice");
  py_assert_ptr(slice);

  gsl_rng* rng= 0;
  
  if(py_random != Py_None) {
    rng= (gsl_rng*) PyCapsule_GetPointer(py_random, "_GSL_RNG");
    py_assert_ptr(rng);
  }

  try {
    cola_lightcones_create(snapshots, sky, remap, slice, lightcones, rng);
  }
  catch (const ColaFileError e) {
    PyErr_SetString(PyExc_IOError, "Unable to read data");
    return NULL;
  }

  Py_RETURN_NONE;
}

