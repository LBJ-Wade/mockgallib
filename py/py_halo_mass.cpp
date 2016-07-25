#include "py_assert.h"
#include "py_halo_mass.h"
#include "halo_mass.h"

static void py_halo_mass_free(PyObject *obj);


PyObject* py_halo_mass_alloc(PyObject* self, PyObject* args)
{
  // _halo_mass_alloc(filename)
  PyObject* bytes;
  char* filename;
  Py_ssize_t len;

  if(!PyArg_ParseTuple(args, "O&", PyUnicode_FSConverter, &bytes)) {
    return NULL;
  }

  PyBytes_AsStringAndSize(bytes, &filename, &len);

  HaloMassFoF* halo_mass;

  try {
    halo_mass= new HaloMassFoF(filename);
  }
  catch(HaloMassFileError) {
    Py_DECREF(bytes);
    PyErr_SetString(PyExc_IOError, "HaloMassFileError");
    return NULL;
  }

  Py_DECREF(bytes);

  return PyCapsule_New(halo_mass, "_HaloMass", py_halo_mass_free);
}


void py_halo_mass_free(PyObject *obj)
{
  HaloMassFoF* const halo_mass=
    (HaloMassFoF*) PyCapsule_GetPointer(obj, "_HaloMass");
  py_assert_void(halo_mass);

  delete halo_mass;
}
