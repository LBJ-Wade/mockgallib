#ifndef PY_HOD_H
#define PY_HOD_H 1

#include "Python.h"

PyObject* py_hod_alloc(PyObject* self, PyObject* args);
PyObject* py_hod_get_coef_all(PyObject* self, PyObject* args);
PyObject* py_hod_get_coef(PyObject* self, PyObject* args);
PyObject* py_hod_set_coef(PyObject* self, PyObject* args);
PyObject* py_hod_ncen(PyObject* self, PyObject* args);
PyObject* py_hod_nsat(PyObject* self, PyObject* args);
PyObject* py_hod_get_z0(PyObject* self, PyObject* args);
PyObject* py_hod_compute_logMmin(PyObject* self, PyObject* args);
PyObject* py_hod_compute_M1(PyObject* self, PyObject* args);
PyObject* py_hod_compute_sigma(PyObject* self, PyObject* args);
PyObject* py_hod_compute_alpha(PyObject* self, PyObject* args);

#endif
