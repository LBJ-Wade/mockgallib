#ifndef PY_MF_H
#define PY_MF_H 1

#include "Python.h"

PyObject* py_mf_alloc(PyObject* self, PyObject* args);
PyObject* py_mf_set_redshift(PyObject* self, PyObject* args);
PyObject* py_mf_f(PyObject* self, PyObject* args);

#endif
