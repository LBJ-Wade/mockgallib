#ifndef PY_CATALOGUE_H
#define PY_CATALOGUE_H 1

#include "Python.h"
#include "lightcone.h"


PyObject* py_lightcones_alloc(PyObject* self, PyObject* args);
PyObject* py_lightcones_load(PyObject* self, PyObject* args);
PyObject* py_lightcones_len(PyObject* self, PyObject* args);

#endif
