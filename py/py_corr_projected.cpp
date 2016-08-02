#include <vector>
#include "corr_projected.h"
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
