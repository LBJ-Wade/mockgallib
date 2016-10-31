#include <iostream>
#include <vector>
#include "corr_projected.h"
#include "hist2d.h"

#include "py_assert.h"
#include "py_corr_projected.h"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"

using namespace std;

PyMODINIT_FUNC
py_corr_projected_module_init()
{
  import_array();

  return NULL;
}

static void py_corr_projected_free(PyObject *obj);
static void py_hist2d_free(PyObject *obj);
  
//
// hist2d
//
PyObject* py_corr_projected_hist2d_alloc(PyObject* self, PyObject* args)
{
  float rp_min, rp_max, pi_max;
  int nbin_rp, nbin_pi;

  if(!PyArg_ParseTuple(args, "ffifi", &rp_min, &rp_max, &nbin_rp,
		       &pi_max, &nbin_pi))
    return NULL;

  Histogram2D<LogBin, LinearBin>* const hist2d=
    new Histogram2D<LogBin, LinearBin>(
	    LogBin(rp_min, rp_max, nbin_rp), LinearBin(0.0f, pi_max, nbin_pi));

  return PyCapsule_New(hist2d, "_Hist2D", py_hist2d_free);  
}

void py_hist2d_free(PyObject *obj)
{
  Histogram2D<LogBin, LinearBin>* const hist2d=
    (Histogram2D<LogBin, LinearBin>*) PyCapsule_GetPointer(obj, "_Hist2D");

  py_assert_void(hist2d);

  delete hist2d;
}

PyObject* py_corr_projected_hist2d_as_array(PyObject* self, PyObject* args)
{
  PyObject *py_hist2d;
  
  if(!PyArg_ParseTuple(args, "O", &py_hist2d)) {
    return NULL;
  }
  
  Histogram2D<LogBin, LinearBin>* const hist2d=
   (Histogram2D<LogBin, LinearBin>*) PyCapsule_GetPointer(py_hist2d, "_Hist2D");
  py_assert_ptr(hist2d);

  const int nd=2;
  npy_intp dims[]= {hist2d->x_nbin(), hist2d->y_nbin()};

  return PyArray_SimpleNewFromData(nd, dims, NPY_DOUBLE, hist2d->hist);
}

PyObject* py_corr_projected_hist2d_set(PyObject* self, PyObject* args)
{
  PyObject *py_hist2d;
  PyObject *bufobj;
  Py_buffer view;
  double npairs;
  
  if(!PyArg_ParseTuple(args, "OOd", &py_hist2d, &bufobj, &npairs))
    return NULL;

  Histogram2D<LogBin, LinearBin>* const hist2d=
   (Histogram2D<LogBin, LinearBin>*) PyCapsule_GetPointer(py_hist2d, "_Hist2D");
  py_assert_ptr(hist2d);

  if(PyObject_GetBuffer(bufobj, &view,
			PyBUF_ANY_CONTIGUOUS | PyBUF_FORMAT) == -1) {
    return NULL;
  }

  if(view.ndim != 2) {
    PyErr_SetString(PyExc_TypeError, "Expected a 2-dimensional array");
    PyBuffer_Release(&view);
    return NULL;
  }

  if(strcmp(view.format,"d") != 0) {
    PyErr_SetString(PyExc_TypeError, "Expected an array of doubles");
    PyBuffer_Release(&view);
    return NULL;
  }

  hist2d->npairs= npairs;

  const int nrow= view.shape[0];
  const int ncol= view.shape[1];

  if(hist2d->x_nbin() != nrow || hist2d->y_nbin() != ncol) {
    PyErr_SetString(PyExc_TypeError, "Unexpected array size");
    PyBuffer_Release(&view);
    return NULL;
  }

  double* p= (double*) view.buf;
  
  const size_t next_row= view.strides[0]/sizeof(double);
  const size_t next_col= view.strides[1]/sizeof(double);

  double* hist= hist2d->hist;
  
  for(int ix=0; ix<nrow; ++ix) {
    for(int iy=0; iy<ncol; ++iy) {  
      hist[ix*ncol + iy]= *(p + next_col*iy);
    }
    p += next_row;
  }

  msg_printf(msg_info, "Set %d rows and %d colums in hist2d; npairs= %e\n", nrow, ncol, hist2d->npairs);

  PyBuffer_Release(&view);
  Py_RETURN_NONE;
}



//
// py_corr_projected
//

PyObject* py_corr_projected_alloc(PyObject* self, PyObject* args)
{
  // _corr_projected_init(rp_min, rp_max, nbin, pi_max, pi_nbin)
  int nbin, pi_nbin;
  float rp_min, rp_max, pi_max;

  if(!PyArg_ParseTuple(args, "ffifi", &rp_min, &rp_max, &nbin,
		       &pi_max, &pi_nbin)) {
    return NULL;
  }
     
  corr_projected_init(rp_min, rp_max, nbin, pi_max, pi_nbin);

  CorrProjected* const corr= new CorrProjected(nbin);

  
  return PyCapsule_New(corr, "_CorrProjected", py_corr_projected_free);  
}


void py_corr_projected_free(PyObject *obj)
{
  CorrProjected* const corr=
    (CorrProjected*) PyCapsule_GetPointer(obj, "_CorrProjected");
  py_assert_void(corr);

  delete corr;
}

PyObject* py_corr_set_radec_min(PyObject* self, PyObject* args)
{
  float ra_min, dec_min;
  
  if(!PyArg_ParseTuple(args, "ff", &ra_min, &dec_min)) {
    return NULL;
  }

  set_radec_min(ra_min, dec_min);
  Py_RETURN_NONE;
}

PyObject* py_corr_projected_compute(PyObject* self, PyObject* args)
{
  // _corr_projected_compute(mock_catalogues, random_catalogues,
  //                         correlation_functions)

  PyObject *py_galaxies, *py_randoms, *py_corr;
  
  if(!PyArg_ParseTuple(args, "OOO", &py_galaxies, &py_randoms, &py_corr)) {
    return NULL;
  }

  Catalogues* const galaxies=
    (Catalogues*) PyCapsule_GetPointer(py_galaxies, "_Catalogues");
  py_assert_ptr(galaxies);

  Catalogues* const randoms=
    (Catalogues*) PyCapsule_GetPointer(py_randoms, "_Catalogues");
  py_assert_ptr(randoms);

  CorrProjected* const corr=
    (CorrProjected*) PyCapsule_GetPointer(py_corr, "_CorrProjected");
  py_assert_ptr(corr);

  corr_projected_compute(galaxies, randoms, corr);

  Py_RETURN_NONE;
}

PyObject* py_corr_projected_compute_rr(PyObject* self, PyObject* args)
{
  // _corr_projected_compute(mock_catalogues, random_catalogues,
  //                         correlation_functions)

  PyObject *py_randoms, *py_hist2d;
  
  if(!PyArg_ParseTuple(args, "OO", &py_randoms, &py_hist2d)) {
    return NULL;
  }
  
  Catalogues* const randoms=
    (Catalogues*) PyCapsule_GetPointer(py_randoms, "_Catalogues");
  py_assert_ptr(randoms);
  

  Histogram2D<LogBin, LinearBin>* const rr=
   (Histogram2D<LogBin, LinearBin>*) PyCapsule_GetPointer(py_hist2d, "_Hist2D");
  py_assert_ptr(rr);

  corr_projected_compute_pairs_rr(randoms, rr);

  return Py_BuildValue("d", rr->npairs);
}

PyObject* py_corr_projected_compute_with_rr(PyObject* self, PyObject* args)
{
  // _corr_projected_compute(mock_catalogues, random_catalogues,
  //                         correlation_functions)

  PyObject *py_galaxies, *py_randoms, *py_hist2d;
  
  if(!PyArg_ParseTuple(args, "OOO", &py_galaxies, &py_randoms, &py_hist2d)) {
    return NULL;
  }

  Catalogues* const galaxies=
    (Catalogues*) PyCapsule_GetPointer(py_galaxies, "_Catalogues");
  py_assert_ptr(galaxies);

  Catalogues* const randoms=
    (Catalogues*) PyCapsule_GetPointer(py_randoms, "_Catalogues");
  py_assert_ptr(randoms);

  Histogram2D<LogBin, LinearBin>* const rr=
   (Histogram2D<LogBin, LinearBin>*) PyCapsule_GetPointer(py_hist2d, "_Hist2D");
  py_assert_ptr(rr);

  corr_projected_compute_with_rr(galaxies, randoms, rr);

  Py_RETURN_NONE;
}


PyObject* py_corr_as_array(PyObject* self, PyObject* args)
{
  PyObject *py_corr;
  
  if(!PyArg_ParseTuple(args, "O", &py_corr)) {
    return NULL;
  }
  CorrProjected* const corr=
    (CorrProjected*) PyCapsule_GetPointer(py_corr, "_CorrProjected");
  py_assert_ptr(corr);

  int nd=2;
  int ncol= 3;
  npy_intp dims[]= {ncol, corr->n};

  return PyArray_SimpleNewFromData(nd, dims, NPY_DOUBLE, corr->rp);
}

PyObject* py_corr_rp_i(PyObject* self, PyObject* args)
{
  // Return rp of vcorr[i] as an array
  int i;
  if(!PyArg_ParseTuple(args, "i", &i)) {
    return NULL;
  }

  CorrProjected* const corr= corr_projected_i(i);

  const int nd=1;
  npy_intp dims[]= {corr->n};

  return PyArray_SimpleNewFromData(nd, dims, NPY_DOUBLE, corr->rp);
}

PyObject* py_corr_wp_i(PyObject* self, PyObject* args)
{
  // Return wp of vcorr[i] as an array  
  int i;
  if(!PyArg_ParseTuple(args, "i", &i)) {
    return NULL;
  }

  CorrProjected* const corr= corr_projected_i(i);

  const int nd=1;
  npy_intp dims[]= {corr->n};

  return PyArray_SimpleNewFromData(nd, dims, NPY_DOUBLE, corr->wp);
}

PyObject* py_corr_rp(PyObject* self, PyObject* args)
{
  PyObject *py_corr;
  
  if(!PyArg_ParseTuple(args, "O", &py_corr)) {
    return NULL;
  }
  CorrProjected* const corr=
    (CorrProjected*) PyCapsule_GetPointer(py_corr, "_CorrProjected");
  py_assert_ptr(corr);

  const int nd=1;
  npy_intp dims[]= {corr->n};

  return PyArray_SimpleNewFromData(nd, dims, NPY_DOUBLE, corr->rp);
}

PyObject* py_corr_wp(PyObject* self, PyObject* args)
{
  PyObject *py_corr;
  
  if(!PyArg_ParseTuple(args, "O", &py_corr)) {
    return NULL;
  }
  CorrProjected* const corr=
    (CorrProjected*) PyCapsule_GetPointer(py_corr, "_CorrProjected");
  py_assert_ptr(corr);

  const int nd=1;
  npy_intp dims[]= {corr->n};

  return PyArray_SimpleNewFromData(nd, dims, NPY_DOUBLE, corr->wp);
}

PyObject* py_corr_dwp(PyObject* self, PyObject* args)
{
  PyObject *py_corr;
  
  if(!PyArg_ParseTuple(args, "O", &py_corr)) {
    return NULL;
  }
  CorrProjected* const corr=
    (CorrProjected*) PyCapsule_GetPointer(py_corr, "_CorrProjected");
  py_assert_ptr(corr);

  const int nd=1;
  npy_intp dims[]= {corr->n};

  return PyArray_SimpleNewFromData(nd, dims, NPY_DOUBLE, corr->dwp);
}
