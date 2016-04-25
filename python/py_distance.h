#ifndef PY_DISTANCE_H
#define PY_DISTANCE_H 1

#include "distance.h"
#include "Python.h"

PyObject* py_distance_init(PyObject* self, PyObject* args);
PyObject* py_distance_redshift(PyObject* self, PyObject* args);

#endif
