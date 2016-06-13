#ifndef PY_NBAR_H
#define PY_NBAR_H 1

#include "Python.h"

PyObject* py_nbar_alloc(PyObject* self, PyObject* args);
PyObject* py_nbar_compute(PyObject* self, PyObject* args);

#endif
