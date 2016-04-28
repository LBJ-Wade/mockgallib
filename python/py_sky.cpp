#include "py_assert.h"
#include "py_sky.h"

static void py_sky_free(PyObject *obj);

PyObject* py_sky_alloc(PyObject* self, PyObject* args)
{
  // _sky_alloc(ra_min, ra_max, dec_min, dec_max, z_min, z_max)
  double ra[2], dec[2], z[2];
  
  if(!PyArg_ParseTuple(args, "dddddd", ra, ra+1, dec, dec+1, z, z+1))
    return NULL;

  Sky* const sky= new Sky(ra, dec, z);

  return PyCapsule_New(sky, "_Sky", py_sky_free);  
}

void py_sky_free(PyObject *obj)
{
  Sky* const sky=
    (Sky*) PyCapsule_GetPointer(obj, "_Sky"); py_assert(sky);

  msg_printf(msg_debug, "freeing sky %x\n", sky);
  delete sky;
}

PyObject* py_sky_box(PyObject* self, PyObject* args)
{
  PyObject* py_sky;
  if(!PyArg_ParseTuple(args, "O", &py_sky))
    return NULL;

  Sky* const sky=
    (Sky*) PyCapsule_GetPointer(py_sky, "_Sky"); py_assert(sky);
  
  
  return Py_BuildValue("(ddd)", sky->width[0], sky->width[1], sky->width[2]);
}
