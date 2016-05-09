#ifndef PY_CATALOGUE_H
#define PY_CATALOGUE_H 1

#include "Python.h"
#include "lightcone.h"

PyMODINIT_FUNC
py_lightcones_module_init();
PyObject* py_lightcones_lighcone(PyObject* self, PyObject* args);

PyObject* py_lightcones_alloc(PyObject* self, PyObject* args);
PyObject* py_lightcones_load_h5(PyObject* self, PyObject* args);
PyObject* py_lightcones_len(PyObject* self, PyObject* args);
PyObject* py_lightcones_lighcone(PyObject* self, PyObject* args);
PyObject* py_lightcones_clear(PyObject* self, PyObject* args);

PyObject* py_lightcone_len(PyObject* self, PyObject* args);
PyObject* py_lightcone_as_array(PyObject* self, PyObject* args);
PyObject* py_lightcone_save_h5(PyObject* self, PyObject* args);
PyObject* py_lightcone_load_h5(PyObject* self, PyObject* args);

#endif
