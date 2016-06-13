#ifndef PY_SIGMA_H
#define PY_SIGMA_H 1

#include "Python.h"

PyMODINIT_FUNC
py_sigma_module_init();

PyObject* py_sigma_alloc(PyObject* self, PyObject* args);
PyObject* py_sigma_n(PyObject* self, PyObject* args);
PyObject* py_sigma_M_range(PyObject* self, PyObject* args);
PyObject* py_sigma_M(PyObject* self, PyObject* args);
PyObject* py_sigma_0inv(PyObject* self, PyObject* args);
PyObject* py_sigma_M_array(PyObject* self, PyObject* args);
PyObject* py_sigma_sinv_array(PyObject* self, PyObject* args);

#endif
