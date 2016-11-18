#ifndef PY_CORR_MULTIPOLE
#define PY_CORR_MULTIPOLE 1

#include "Python.h"

PyMODINIT_FUNC
py_corr_multipole_module_init();

PyObject* py_corr_multipole_hist2d_alloc(PyObject* self, PyObject* args);
PyObject* py_corr_multipole_hist2d_as_array(PyObject* self, PyObject* args);
PyObject* py_corr_multipole_hist2d_set(PyObject* self, PyObject* args);

PyObject* py_corr_multipole_alloc(PyObject* self, PyObject* args);
PyObject* py_corr_multipole_set_radec_min(PyObject* self, PyObject* args);
PyObject* py_corr_multipole_set_pair_correction(PyObject* self, PyObject* args);

PyObject* py_corr_multipole_compute(PyObject* self, PyObject* args);
PyObject* py_corr_multipole_as_array(PyObject* self, PyObject* args);

PyObject* py_corr_multipole_compute_rr(PyObject* self, PyObject* args);
PyObject* py_corr_multipole_compute_with_rr(PyObject* self, PyObject* args);


PyObject* py_corr_multipole_r_i(PyObject* self, PyObject* args);
PyObject* py_corr_multipole_xi0_i(PyObject* self, PyObject* args);
PyObject* py_corr_multipole_xi2_i(PyObject* self, PyObject* args);
#endif
