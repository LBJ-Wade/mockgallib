#ifndef PY_COMM_H
#define PY_COMM_H 1

#include "Python.h"

PyObject* py_comm_init(PyObject *self, PyObject* args);
PyObject* py_comm_finalise(PyObject *self, PyObject* args);
PyObject* py_comm_this_rank(PyObject *self, PyObject* args);
PyObject* py_comm_n_nodes(PyObject *self, PyObject* args);
PyObject* py_comm_bcast_str(PyObject *self, PyObject* args);
#endif
