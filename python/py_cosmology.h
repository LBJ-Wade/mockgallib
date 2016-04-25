#ifndef PY_COSMOLOGY_H
#define PY_COSMOLOGY_H 1

#include "Python.h"

PyObject* py_cosmology_set(PyObject* self, PyObject* args);
PyObject* py_cosmology_rhom(PyObject* self, PyObject* args);
PyObject* py_cosmology_compute_comoving_distance(PyObject* self, PyObject* args);
#endif
