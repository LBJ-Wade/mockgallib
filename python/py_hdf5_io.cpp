#include "lightcone.h"
#include "hdf5_io.h"
#include "py_assert.h"
#include "py_hdf5_io.h"

PyObject* py_hdf5_write_lightcone(PyObject* self, PyObject* args)
{
  // _hdf5_write_lightcone(filename, _lightcone)
  PyObject* py_lightcone;
  PyObject* bytes;
  char* filename;
  Py_ssize_t len;

  if(!PyArg_ParseTuple(args, "O&O", PyUnicode_FSConverter, &bytes,
		       &py_lightcone)) {
    return NULL;
  }

  PyBytes_AsStringAndSize(bytes, &filename, &len);

  LightCone* const lightcone=
    (LightCone *) PyCapsule_GetPointer(py_lightcone, "_LightCone");
  py_assert_ptr(lightcone);

  hdf5_write_lightcone(filename, lightcone);

  Py_DECREF(bytes);
  Py_RETURN_NONE;
}

