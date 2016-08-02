#ifndef PY_CORR_PROJECTED
#define PY_CORR_PROJECTED 1

#include "Python.h"

PyMODINIT_FUNC
py_corr_projected_module_init();

PyObject* py_corr_projected_alloc(PyObject* self, PyObject* args);
PyObject* py_corr_projected_compute(PyObject* self, PyObject* args);
PyObject* py_corr_as_array(PyObject* self, PyObject* args);

PyObject* py_corr_rp(PyObject* self, PyObject* args);
PyObject* py_corr_wp(PyObject* self, PyObject* args);
PyObject* py_corr_dwp(PyObject* self, PyObject* args);
#endif
