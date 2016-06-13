#ifndef PY_SKY_H
#define PY_SKY_H 1

#include "Python.h"
#include "sky.h"

PyObject* py_sky_alloc(PyObject* self, PyObject* args);
PyObject* py_sky_boxsize(PyObject* self, PyObject* args);
PyObject* py_sky_left(PyObject* self, PyObject* args);
PyObject* py_sky_right(PyObject* self, PyObject* args);
PyObject* py_sky_r_range(PyObject* self, PyObject* args);

PyObject* py_sky_compute_radec(PyObject* self, PyObject* args);
PyObject* py_sky_compute_x(PyObject* self, PyObject* args);
#endif
