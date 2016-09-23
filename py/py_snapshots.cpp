#include <stdexcept>
#include "py_assert.h"
#include "snapshot.h"
#include "py_snapshots.h"

static void py_snapshots_free(PyObject *obj);


PyObject* py_snapshots_alloc(PyObject* self, PyObject* args)
{
  Snapshots* snps= new Snapshots();

  return PyCapsule_New(snps, "_Snapshots", py_snapshots_free);
}

void py_snapshots_free(PyObject *obj)
{
  Snapshots* const snps=
    (Snapshots*) PyCapsule_GetPointer(obj, "_Snapshots");
  py_assert_void(snps);

  delete snps;
}

PyObject* py_snapshots_insert(PyObject* self, PyObject* args)
{
  // _snapshots_insert(_snps, fof_filename, halo_mass_filename,
  //                   M_part_min, M_halo_min,
  //                   a_snp, a_min, a_max, _halo_mass)
  PyObject *py_snps, *bytes_fof, *bytes_part, *bytes_halo_mass;
  double a_snp, a_min, a_max;
  double M_part_min, M_halo_min;
  
  if(!PyArg_ParseTuple(args, "OO&O&O&ddddd", &py_snps,
		       PyUnicode_FSConverter, &bytes_fof,
		       PyUnicode_FSConverter, &bytes_part,
		       PyUnicode_FSConverter, &bytes_halo_mass,
		       &M_part_min, &M_halo_min,
		       &a_snp, &a_min, &a_max)) {
    return NULL;
  }

  char *filename_fof, *filename_part, *filename_halo_mass;
  Py_ssize_t len;

  PyBytes_AsStringAndSize(bytes_fof, &filename_fof, &len);
  PyBytes_AsStringAndSize(bytes_part, &filename_part, &len);
  PyBytes_AsStringAndSize(bytes_halo_mass, &filename_halo_mass, &len);

  Snapshots* snps= (Snapshots*) PyCapsule_GetPointer(py_snps, "_Snapshots");
  py_assert_ptr(snps);

  Snapshot* snp;
  try {
    snp= new Snapshot(filename_fof, filename_part, filename_halo_mass,
		      M_part_min, M_halo_min,
		      a_snp, a_min, a_max);
  }
  catch(const HaloMassFileError) {
    PyErr_SetString(PyExc_IOError, "Unable to read halo mass file");
    return NULL;
  }

  snps->push_back(snp);
  

  Py_DECREF(bytes_fof);
  Py_DECREF(bytes_part);
  Py_DECREF(bytes_halo_mass);

  Py_RETURN_NONE;
}

PyObject* py_snapshots_len(PyObject* self, PyObject* args)
{
  // _snapshots_get(_snps)
  PyObject *py_snps;
  
  if(!PyArg_ParseTuple(args, "O", &py_snps)) {
    return NULL;
  }

  Snapshots* snps= (Snapshots*) PyCapsule_GetPointer(py_snps, "_Snapshots");
  py_assert_ptr(snps);

  return Py_BuildValue("i", (int) snps->size());
}

PyObject* py_snapshots_get(PyObject* self, PyObject* args)
{
  // _snapshots_get(_snps, i)
  PyObject *py_snps;
  int i;
  
  if(!PyArg_ParseTuple(args, "Oi", &py_snps, &i)) {
    return NULL;
  }

  Snapshots* snps= (Snapshots*) PyCapsule_GetPointer(py_snps, "_Snapshots");
  py_assert_ptr(snps);

  Snapshot const* snp;
  try {
    snp= snps->at(i);
  }
  catch(const std::out_of_range) {
    PyErr_SetNone(PyExc_IndexError);
    return NULL;
  }

  return Py_BuildValue("(sddd)", snp->filename_fof,
		       snp->a_snp, snp->a_min, snp->a_max);
}

PyObject* py_snapshots_clear(PyObject* self, PyObject* args)
{
  // _snapshots_get(_snps)
  PyObject *py_snps;
  
  if(!PyArg_ParseTuple(args, "O", &py_snps)) {
    return NULL;
  }

  Snapshots* snps= (Snapshots*) PyCapsule_GetPointer(py_snps, "_Snapshots");
  py_assert_ptr(snps);

  snps->clear();

  Py_RETURN_NONE;
}
