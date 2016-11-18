#include <iostream>
#include <vector>
#include "corr_multipole.h"
#include "hist2d.h"

#include "py_assert.h"
#include "py_corr_multipole.h"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"

using namespace std;

PyMODINIT_FUNC
py_corr_multipole_module_init()
{
  import_array();

  return NULL;
}

static void py_corr_multipole_free(PyObject *obj);
static void py_hist2d_free(PyObject *obj);
  
//
// hist2d
//
PyObject* py_corr_multipole_hist2d_alloc(PyObject* self, PyObject* args)
{
  float r_min, r_max;
  int nbin_r, nbin_mu;

  if(!PyArg_ParseTuple(args, "ffii", &r_min, &r_max, &nbin_r, &nbin_mu))
    return NULL;

  Histogram2D<LogBin, LinearBin>* const hist2d=
    new Histogram2D<LogBin, LinearBin>(
	    LogBin(r_min, r_max, nbin_r), LinearBin(0.0f, 1.0, nbin_mu));

  return PyCapsule_New(hist2d, "_Hist2D", py_hist2d_free);  
}

void py_hist2d_free(PyObject *obj)
{
  Histogram2D<LogBin, LinearBin>* const hist2d=
    (Histogram2D<LogBin, LinearBin>*) PyCapsule_GetPointer(obj, "_Hist2D");

  py_assert_void(hist2d);

  delete hist2d;
}

PyObject* py_corr_multipole_hist2d_as_array(PyObject* self, PyObject* args)
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

PyObject* py_corr_multipole_hist2d_set(PyObject* self, PyObject* args)
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
// py_corr_multipole
//

PyObject* py_corr_multipole_alloc(PyObject* self, PyObject* args)
{
  // _corr_multipole_init(r_min, r_max, nbin, nbin_mu)
  int nbin_r, nbin_mu;
  float r_min, r_max;

  if(!PyArg_ParseTuple(args, "ffii", &r_min, &r_max, &nbin_r, &nbin_mu)) {
    return NULL;
  }
     
  corr_multipole_init(r_min, r_max, nbin_r, nbin_mu);

  CorrMultipole* const corr= new CorrMultipole(nbin_r);

  
  return PyCapsule_New(corr, "_CorrMultipole", py_corr_multipole_free);  
}


void py_corr_multipole_free(PyObject *obj)
{
  CorrMultipole* const corr=
    (CorrMultipole*) PyCapsule_GetPointer(obj, "_CorrMultipole");
  py_assert_void(corr);

  delete corr;
}

PyObject* py_corr_multipole_set_radec_min(PyObject* self, PyObject* args)
{
  float ra_min, dec_min;
  
  if(!PyArg_ParseTuple(args, "ff", &ra_min, &dec_min)) {
    return NULL;
  }

  corr_multipole_set_radec_min(ra_min, dec_min);
  Py_RETURN_NONE;
}

PyObject* py_corr_multipole_set_pair_correction(PyObject* self, PyObject* args)
{
  PyObject* bytes;
  
  if(!PyArg_ParseTuple(args, "O&", PyUnicode_FSConverter, &bytes)) {
    return NULL;
  }
  
  char* filename;
  Py_ssize_t len;
  PyBytes_AsStringAndSize(bytes, &filename, &len);

  corr_multipole_set_pair_correction(filename);

  Py_DECREF(bytes);
  
  Py_RETURN_NONE;
}

  

PyObject* py_corr_multipole_compute(PyObject* self, PyObject* args)
{
  // _corr_multipole_compute(mock_catalogues, random_catalogues,
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

  CorrMultipole* const corr=
    (CorrMultipole*) PyCapsule_GetPointer(py_corr, "_CorrMultipole");
  py_assert_ptr(corr);

  corr_multipole_compute(galaxies, randoms);

  Py_RETURN_NONE;
}

PyObject* py_corr_multipole_compute_rr(PyObject* self, PyObject* args)
{
  // _corr_multipole_compute(mock_catalogues, random_catalogues,
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

  corr_multipole_compute_pairs_rr(randoms, rr);

  return Py_BuildValue("d", rr->npairs);
}

PyObject* py_corr_multipole_compute_with_rr(PyObject* self, PyObject* args)
{
  // _corr_multipole_compute(mock_catalogues, random_catalogues,
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

  corr_multipole_compute_with_rr(galaxies, randoms, rr);

  Py_RETURN_NONE;
}


PyObject* py_corr_multipole_as_array(PyObject* self, PyObject* args)
{
  PyObject *py_corr;
  
  if(!PyArg_ParseTuple(args, "O", &py_corr)) {
    return NULL;
  }
  CorrMultipole* const corr=
    (CorrMultipole*) PyCapsule_GetPointer(py_corr, "_CorrMultipole");
  py_assert_ptr(corr);

  int nd=2;
  int ncol= 3;
  npy_intp dims[]= {ncol, corr->n};

  return PyArray_SimpleNewFromData(nd, dims, NPY_DOUBLE, corr->r);
}

PyObject* py_corr_multipole_r_i(PyObject* self, PyObject* args)
{
  // Return rp of vcorr[i] as an array
  int i;
  if(!PyArg_ParseTuple(args, "i", &i)) {
    return NULL;
  }

  CorrMultipole* const corr= corr_multipole_i(i);

  const int nd=1;
  npy_intp dims[]= {corr->n};

  return PyArray_SimpleNewFromData(nd, dims, NPY_DOUBLE, corr->r);
}

PyObject* py_corr_multipole_xi0_i(PyObject* self, PyObject* args)
{
  // Return wp of vcorr[i] as an array  
  int i;
  if(!PyArg_ParseTuple(args, "i", &i)) {
    return NULL;
  }

  CorrMultipole* const corr= corr_multipole_i(i);

  const int nd=1;
  npy_intp dims[]= {corr->n};

  return PyArray_SimpleNewFromData(nd, dims, NPY_DOUBLE, corr->xi0);
}

PyObject* py_corr_multipole_xi2_i(PyObject* self, PyObject* args)
{
  // Return wp of vcorr[i] as an array  
  int i;
  if(!PyArg_ParseTuple(args, "i", &i)) {
    return NULL;
  }

  CorrMultipole* const corr= corr_multipole_i(i);

  const int nd=1;
  npy_intp dims[]= {corr->n};

  return PyArray_SimpleNewFromData(nd, dims, NPY_DOUBLE, corr->xi2);
}

