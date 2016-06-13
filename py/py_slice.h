#ifndef PY_SLICE_H
#define PY_SLICE_H 1

#include "Python.h"

PyObject* py_slice_alloc(PyObject* self, PyObject* args);
PyObject* py_slice_len(PyObject* self, PyObject* args);
PyObject* py_slice_boxsize(PyObject* self, PyObject* args);

#endif
