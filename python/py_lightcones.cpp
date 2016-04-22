#include "msg.h"
#include "lightcone.h"
#include "py_lightcones.h"
#include "hdf5_io.h"

using namespace std;

static void py_lightcones_free(PyObject *obj);

PyObject* py_lightcones_alloc(PyObject* self, PyObject* args)
{
  LightCones* const lightcones=
    new LightCones();

  return PyCapsule_New(lightcones, "_LightCones", py_lightcones_free);  
}

void py_lightcones_free(PyObject *obj)
{
  LightCones* const lightcones=
    (LightCones*) PyCapsule_GetPointer(obj, "_LightCones");

  msg_printf(msg_debug, "freeing lightcones %x\n", lightcones);

  delete lightcones;
}

PyObject* py_lightcones_load(PyObject* self, PyObject* args)
{
  PyObject* py_lightcones;
  PyObject* bytes;
  char* filename;
  Py_ssize_t len;

  if(!PyArg_ParseTuple(args, "OO&", &py_lightcones,
		       PyUnicode_FSConverter, &bytes)) {
    return NULL;
  }
  
  PyBytes_AsStringAndSize(bytes, &filename, &len);

  LightCone* const lightcone= new LightCone();

  try {
    hdf5_read_lightcone(filename, lightcone);
  }
  catch(LightconeFileError) {
    delete lightcone;
    Py_DECREF(bytes);
    PyErr_SetString(PyExc_IOError, "LightconeFileError");
    Py_RETURN_NONE;
  }

  Py_DECREF(bytes);

  LightCones* const lightcones=
    (LightCones*) PyCapsule_GetPointer(py_lightcones, "_LightCones");

  lightcones->push_back(lightcone);

  Py_RETURN_NONE;
}

PyObject* py_lightcones_len(PyObject* self, PyObject* args)
{
  PyObject* py_lightcones;
  if(!PyArg_ParseTuple(args, "O", &py_lightcones)) {
    return NULL;
  }

  LightCones* const lightcones=
    (LightCones*) PyCapsule_GetPointer(py_lightcones, "_LightCones");

  return Py_BuildValue("i", (int) lightcones->size());
}
