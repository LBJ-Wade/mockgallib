#ifndef PY_POWER_H
#define PY_POWER_H 1

#include <vector>
#include "Python.h"

PyMODINIT_FUNC
py_power_module_init();
//void py_power_add_methods(std::vector<PyMethodDef>& methods);
PyObject* py_power_alloc(PyObject* self, PyObject* args);
PyObject* py_power_sigma(PyObject* self, PyObject* args);
PyObject* py_power_n(PyObject* self, PyObject* args);
PyObject* py_power_k(PyObject* self, PyObject* args);
PyObject* py_power_P(PyObject* self, PyObject* args);
PyObject* py_power_ki(PyObject* self, PyObject* args);
PyObject* py_power_Pi(PyObject* self, PyObject* args);

#endif
