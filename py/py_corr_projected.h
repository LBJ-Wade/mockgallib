#ifndef PY_CORR_PROJECTED
#define PY_CORR_PROJECTED 1

#include "Python.h"

PyMODINIT_FUNC
py_corr_projected_module_init();

PyObject* py_corr_projected_hist2d_alloc(PyObject* self, PyObject* args);
PyObject* py_corr_projected_hist2d_as_array(PyObject* self, PyObject* args);
PyObject* py_corr_projected_hist2d_set(PyObject* self, PyObject* args);

PyObject* py_corr_projected_alloc(PyObject* self, PyObject* args);
PyObject* py_corr_set_radec_min(PyObject* self, PyObject* args);
PyObject* py_corr_set_pair_correction(PyObject* self, PyObject* args);

PyObject* py_corr_projected_compute(PyObject* self, PyObject* args);
PyObject* py_corr_as_array(PyObject* self, PyObject* args);

PyObject* py_corr_rp(PyObject* self, PyObject* args);
PyObject* py_corr_wp(PyObject* self, PyObject* args);
PyObject* py_corr_dwp(PyObject* self, PyObject* args);

PyObject* py_corr_rp_i(PyObject* self, PyObject* args);
PyObject* py_corr_wp_i(PyObject* self, PyObject* args);
PyObject* py_corr_projected_compute_rr(PyObject* self, PyObject* args);
PyObject* py_corr_projected_compute_with_rr(PyObject* self, PyObject* args);
#endif
