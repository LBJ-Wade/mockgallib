#ifndef PY_CATALOGUES_H
#define PY_CATALOGUES_H 1

#include "Python.h"

PyMODINIT_FUNC
py_catalogues_module_init();

PyObject* py_catalogue_len(PyObject* self, PyObject* args);

PyObject* py_catalogues_alloc(PyObject* self, PyObject* args);
PyObject* py_catalogues_load_h5(PyObject* self, PyObject* args);

PyObject* py_catalogues_generate(PyObject* self, PyObject* args);
PyObject* py_catalogues_len(PyObject* self, PyObject* args);
PyObject* py_catalogues_catalogue(PyObject* self, PyObject* args);
PyObject* py_catalogues_append(PyObject* self, PyObject* args);

PyObject* py_catalogues_ngal(PyObject* self, PyObject* args);
PyObject* py_catalogues_nz(PyObject* self, PyObject* args);
PyObject* py_catalogues_subsample(PyObject* self, PyObject* args);
#endif

