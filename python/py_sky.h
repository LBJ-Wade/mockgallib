#ifndef PY_SKY_H
#define PY_SKY_H 1

#include "Python.h"
#include "sky.h"

PyObject* py_sky_alloc(PyObject* self, PyObject* args);
PyObject* py_sky_box(PyObject* self, PyObject* args);

#endif
