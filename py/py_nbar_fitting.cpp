#include "msg.h"
#include "nbar_fitting.h"
#include "py_nbar_fitting.h"
#include <iostream>

using namespace std;


static void py_nbar_fitting_free(PyObject *obj);

PyObject* py_nbar_fitting_alloc(PyObject* self, PyObject* args)
{
  // _nbar_fitting_alloc(_ps, _hod, _nbar_obs, z_min, z_max)

  double z_min, z_max;
  PyObject *py_hod;
  PyObject *bufobj;
  Py_buffer view;

  
  if(!PyArg_ParseTuple(args, "OOdd", &py_hod, &bufobj, &z_min, &z_max))
    return NULL;

  Hod* const hod=
    (Hod*) PyCapsule_GetPointer(py_hod, "_HOD");

  if(PyObject_GetBuffer(bufobj, &view,
			PyBUF_ANY_CONTIGUOUS | PyBUF_FORMAT) == -1) {
    return NULL;
  }

  //
  // Read array and copy to vector<Nbar>
  //
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

  const int n= view.shape[0];
  const int ncol= view.shape[1];

  if(n == 0) {
    PyErr_SetString(PyExc_TypeError, "no data in nbar_obs array");
    return NULL;
  }
  
  if(ncol < 2) {
    PyErr_SetString(PyExc_TypeError, "Expected two columns for z and nbar");
    return NULL;
  }
  else if (ncol > 3) {
    msg_printf(msg_warn, "Neglecting 4th and later columns in nbar_obs data for nbar_fitting\n");
  }

  double* p= (double*) view.buf;
  
  const size_t next_row= view.strides[0]/sizeof(double);
  const size_t next_col= view.strides[1]/sizeof(double);

  Nbar nbar;
  vector<Nbar>* v= new vector<Nbar>;

  if(ncol >= 3) {
    for(int i=0; i<n; ++i) {
      nbar.z= *p;
      nbar.nbar=  *(p + next_col);
      nbar.dnbar= *(p + 2*next_col);
      
      if(z_min <= nbar.z && nbar.z <= z_max) {
	v->push_back(nbar);
	//printf("%e %e %e\n", nbar.z, nbar.nbar, nbar.dnbar);
      }
      p += next_row;
    }
    msg_printf(msg_info, "Using %lu data of z, nbar, dnbar for nbar fitting\n", v->size()); 
  }
  else {
    for(int i=0; i<n; ++i) {
      nbar.z= *p;
      nbar.nbar= *(p + next_col);
      nbar.dnbar= nbar.nbar; // minimising the relative error
      
      if(z_min <= nbar.z && nbar.z <= z_max) {
	v->push_back(nbar);
	//printf("%e %e\n", nbar.z, nbar.nbar);
      }
      p += next_row;
    }
    msg_printf(msg_info, "Using %lu data of z nbar for nbar fitting\n", v->size()); 

  }

  
  NbarFitting* const fitting=
    nbar_fitting_alloc(hod, v, z_min, z_max);

  return PyCapsule_New(fitting, "_NbarFitting", py_nbar_fitting_free);
}


void py_nbar_fitting_free(PyObject *obj)
{
  NbarFitting* const fitting=
      (NbarFitting*) PyCapsule_GetPointer(obj, "_NbarFitting");
  msg_printf(msg_debug, "Freeing nbar_fitting\n");

  nbar_fitting_free(fitting);
}

PyObject* py_nbar_fitting_len(PyObject* self, PyObject* args)
{
  PyObject* py_fitting;
  if(!PyArg_ParseTuple(args, "O", &py_fitting))
    return NULL;

  NbarFitting* const fitting=
    (NbarFitting*) PyCapsule_GetPointer(py_fitting, "_NbarFitting");

  return Py_BuildValue("i", fitting->vobs->size());
}

PyObject* py_nbar_fitting_z(PyObject* self, PyObject* args)
{
  // return vobs[i].z
  int i;
  PyObject* py_fitting;
  if(!PyArg_ParseTuple(args, "Oi", &py_fitting, &i))
    return NULL;

  NbarFitting* const fitting=
    (NbarFitting*) PyCapsule_GetPointer(py_fitting, "_NbarFitting");

  //cerr << i << " " << fitting->vobs->at(i).z << endl;

  return Py_BuildValue("d", fitting->vobs->at(i).z);
}

PyObject* py_nbar_fitting_nbar_obs(PyObject* self, PyObject* args)
{
  // return vobs[i].nbar
  int i;
  PyObject* py_fitting;
  if(!PyArg_ParseTuple(args, "Oi", &py_fitting, &i))
    return NULL;

  NbarFitting* const fitting=
    (NbarFitting*) PyCapsule_GetPointer(py_fitting, "_NbarFitting");

  return Py_BuildValue("d", fitting->vobs->at(i).nbar);
}

PyObject* py_nbar_fitting_nbar_hod(PyObject* self, PyObject* args)
{
  // return vhod[i].nbar
  int i;
  PyObject* py_fitting;
  if(!PyArg_ParseTuple(args, "Oi", &py_fitting, &i))
    return NULL;

  NbarFitting* const fitting=
    (NbarFitting*) PyCapsule_GetPointer(py_fitting, "_NbarFitting");

  return Py_BuildValue("d", fitting->vhod->at(i).nbar);
}

PyObject* py_nbar_fitting_hod(PyObject* self, PyObject* args)
{
  PyObject* py_fitting;
  if(!PyArg_ParseTuple(args, "O", &py_fitting))
    return NULL;

  NbarFitting* const fitting=
    (NbarFitting*) PyCapsule_GetPointer(py_fitting, "_NbarFitting");

  //printf("%e %e %e\n", fitting->hod->c[0], fitting->hod->c[1], fitting->hod->c[2]);

  
  return PyCapsule_New(fitting->hod, "_HOD", NULL);
}

PyObject* py_nbar_fitting_compute(PyObject* self, PyObject* args)
{
  PyObject* py_fitting;
  if(!PyArg_ParseTuple(args, "O", &py_fitting))
    return NULL;

  NbarFitting* const fitting=
    (NbarFitting*) PyCapsule_GetPointer(py_fitting, "_NbarFitting");

  nbar_fitting_compute(fitting);

  cerr << "fitting_compute " << fitting->iter << " " << fitting->chi2 << endl;

  return Py_BuildValue("(id)", fitting->iter, fitting->chi2);
}
