#ifndef PY_CALLBACK_H
#define PY_CALLBACK_H 1

#include "Python.h"

PyObject* py_callback_standby(PyObject *self, PyObject *args);
PyObject* py_callback_sync(PyObject *self, PyObject *args);
PyObject* py_callback_release(PyObject *self, PyObject *args);

#endif
