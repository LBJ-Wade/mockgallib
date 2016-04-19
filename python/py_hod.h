#ifndef PY_HOD_H
#define PY_HOD_H 1

#include "Python.h"

PyObject* py_hod_alloc(PyObject* self, PyObject* args);
PyObject* py_hod_get_coef(PyObject* self, PyObject* args);
PyObject* py_hod_ncen(PyObject* self, PyObject* args);
PyObject* py_hod_nsat(PyObject* self, PyObject* args);

#endif
