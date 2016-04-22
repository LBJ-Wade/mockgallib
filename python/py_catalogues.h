#ifndef PY_CATALOGUES_H
#define PY_CATALOGUES_H 1

#include "Python.h"
#include "catalogue.h"

PyObject* py_catalogues_alloc(PyObject* self, PyObject* args);
PyObject* py_catalogues_generate(PyObject* self, PyObject* args);
PyObject* py_catalogues_len(PyObject* self, PyObject* args);

#endif

