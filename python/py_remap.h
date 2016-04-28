#ifndef PY_REMAP_H
#define PY_REMAP_H 1

#include "remap.h"
#include "Python.h"

PyObject* py_remap_alloc(PyObject* self, PyObject* args);
PyObject* py_remap_boxsize(PyObject* self, PyObject* args);

#endif
