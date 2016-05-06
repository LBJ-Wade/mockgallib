#include <iostream>
#include "py_assert.h"
#include "py_sky.h"

using namespace std; // debug !!!

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
    (Sky*) PyCapsule_GetPointer(obj, "_Sky"); py_assert_void(sky);

  msg_printf(msg_debug, "freeing sky %x\n", sky);
  delete sky;
}

PyObject* py_sky_boxsize(PyObject* self, PyObject* args)
{
  PyObject* py_sky;
  if(!PyArg_ParseTuple(args, "O", &py_sky))
    return NULL;

  Sky* const sky=
    (Sky*) PyCapsule_GetPointer(py_sky, "_Sky"); py_assert_ptr(sky);
  
  
  return Py_BuildValue("(ddd)", sky->width[0], sky->width[1], sky->width[2]);
}

PyObject* py_sky_left(PyObject* self, PyObject* args)
{
  PyObject* py_sky;
  if(!PyArg_ParseTuple(args, "O", &py_sky))
    return NULL;

  Sky* const sky=
    (Sky*) PyCapsule_GetPointer(py_sky, "_Sky"); py_assert_ptr(sky);
  
  
  return Py_BuildValue("(ddd)", sky->left[0], sky->left[1], sky->left[2]);
}

PyObject* py_sky_right(PyObject* self, PyObject* args)
{
  PyObject* py_sky;
  if(!PyArg_ParseTuple(args, "O", &py_sky))
    return NULL;

  Sky* const sky=
    (Sky*) PyCapsule_GetPointer(py_sky, "_Sky"); py_assert_ptr(sky);

  
  
  
  return Py_BuildValue("(ddd)", sky->right[0], sky->right[1], sky->right[2]);
}

PyObject* py_sky_r_range(PyObject* self, PyObject* args)
{
  PyObject* py_sky;
  if(!PyArg_ParseTuple(args, "O", &py_sky))
    return NULL;

  Sky* const sky=
    (Sky*) PyCapsule_GetPointer(py_sky, "_Sky"); py_assert_ptr(sky);

  return Py_BuildValue("(dd)", sky->r_range[0], sky->r_range[1]);
}

PyObject* py_sky_compute_radec(PyObject* self, PyObject* args)
{
  // _sky_compute_radec(_sky, x, y, z)
  PyObject* py_sky;
  float x[3];
  if(!PyArg_ParseTuple(args, "Offf", &py_sky, x, x+1, x+2))
    return NULL;

  Sky* const sky=
    (Sky*) PyCapsule_GetPointer(py_sky, "_Sky"); py_assert_ptr(sky);

  float radec[2];
  sky->compute_radec(x, radec);
  
  return Py_BuildValue("(dd)", (double)radec[0], (double)radec[1]);
}

PyObject* py_sky_compute_x(PyObject* self, PyObject* args)
{
  // _sky_compute_x(_sky, r, ra, dec)
  PyObject* py_sky;
  float r, radec[2];
  if(!PyArg_ParseTuple(args, "Offf", &py_sky, &r, radec, radec+1))
    return NULL;

  Sky* const sky=
    (Sky*) PyCapsule_GetPointer(py_sky, "_Sky"); py_assert_ptr(sky);

  float x[3];
  sky->compute_x(r, radec, x);
  
  return Py_BuildValue("(ddd)", (double) x[0], (double) x[1], (double) x[2]);
}

