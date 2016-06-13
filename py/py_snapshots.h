#include "Python.h"
#include "snapshot.h"

PyObject* py_snapshots_alloc(PyObject* self, PyObject* args);
PyObject* py_snapshots_insert(PyObject* self, PyObject* args);
PyObject* py_snapshots_len(PyObject* self, PyObject* args);
PyObject* py_snapshots_get(PyObject* self, PyObject* args);
PyObject* py_snapshots_clear(PyObject* self, PyObject* args);
