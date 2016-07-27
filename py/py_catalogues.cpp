#include <iostream>
#include <stdexcept>
#include "catalogue.h"
#include "py_assert.h"
#include "py_catalogues.h"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"

using namespace std;

PyMODINIT_FUNC
py_catalogues_module_init()
{
  import_array();

  return NULL;
}

static void py_catalogues_free(PyObject *obj);


PyObject* py_catalogues_alloc(PyObject* self, PyObject* args)
{
  Catalogues* const cats= new Catalogues();

  return PyCapsule_New(cats, "_Catalogues", py_catalogues_free);  
}

void py_catalogues_free(PyObject *obj)
{
  Catalogues* const cats=
    (Catalogues*) PyCapsule_GetPointer(obj, "_Catalogues");
  py_assert_void(cats);

  msg_printf(msg_debug, "freeing catalogues %x\n", cats);
  delete cats;
}

PyObject* py_catalogues_generate(PyObject* self, PyObject* args)
{
  // _catalogues_generate_galaxies(_catalogues, _hod, _lightcone, z_min, z_max, random)
  // random: 0 for mock catalogue, 1 for randoms
  
  PyObject *py_catalogues, *py_hod, *py_lightcones;
  double z_min, z_max;
  int random;
  if(!PyArg_ParseTuple(args, "OOOddi", &py_catalogues, &py_hod, &py_lightcones,
		       &z_min, &z_max, &random)) {
    return NULL;
  }

  Catalogues* const cats=
    (Catalogues*) PyCapsule_GetPointer(py_catalogues, "_Catalogues");
  py_assert_ptr(cats);

  Hod* const hod=
    (Hod*) PyCapsule_GetPointer(py_hod, "_HOD");
  py_assert_ptr(hod);

  LightCones* const lightcones=
    (LightCones*) PyCapsule_GetPointer(py_lightcones, "_LightCones");
  py_assert_ptr(lightcones);

  if(cats->empty())
    cats->allocate(lightcones->size());

  py_assert_ptr(cats->size() == lightcones->size());


  const size_t n= lightcones->size();

  if(random) {
    for(size_t i=0; i<n; ++i) {
      catalogue_generate_random(hod, lightcones->at(i),
				z_min, z_max, cats->at(i));
      msg_printf(msg_verbose,
		 "random catalogue %lu generated with %lu galaxies\n",
		 i, cats->at(i)->size());
      
    }
  }
  else {
    for(size_t i=0; i<n; ++i) {
      catalogue_generate_mock(hod, lightcones->at(i),
			      z_min, z_max, cats->at(i));
      
      msg_printf(msg_verbose,
		 "mock catalogue %lu generated with %lu galaxies\n",
		 i, cats->at(i)->size());
    }
  }

  Py_RETURN_NONE;
}


PyObject* py_catalogues_len(PyObject* self, PyObject* args)
{
  // return number of catalogue in the catalogues
  PyObject* py_catalogues;
  if(!PyArg_ParseTuple(args, "O", &py_catalogues)) {
    return NULL;
  }

  Catalogues* const catalogues=
    (Catalogues*) PyCapsule_GetPointer(py_catalogues, "_Catalogues");
  py_assert_ptr(catalogues);

  return Py_BuildValue("i", (int) catalogues->size());
}


PyObject* py_catalogues_catalogue(PyObject* self, PyObject* args)
{
  // _catalogues_catalogue(_cats, i)
  // return ith catalogue as an np.array
  
  PyObject* py_catalogues;
  int i;
  if(!PyArg_ParseTuple(args, "Oi", &py_catalogues, &i)) {
    return NULL;
  }

  Catalogues* const catalogues=
    (Catalogues*) PyCapsule_GetPointer(py_catalogues, "_Catalogues");
  py_assert_ptr(catalogues);

  py_assert_ptr(sizeof(Particle) % sizeof(double) == 0);

  Catalogue* cat= 0;
  try {
    cat= catalogues->at(i);
  }
  catch(const std::out_of_range) {
    PyErr_SetNone(PyExc_IndexError);
    return NULL;
  }

  
  int nd=2;
  int ncol= sizeof(Particle)/sizeof(float);
  npy_intp dims[]= {(npy_intp) cat->size(), ncol};

  return PyArray_SimpleNewFromData(nd, dims, NPY_FLOAT, &(cat->front()));
}


PyObject* py_catalogues_append(PyObject* self, PyObject* args)
{
  // _catalogues_append(_catalogues, array)
  // array must have 4 columns, x,y,z,w
  
  PyObject *py_catalogues, *py_array;
  if(!PyArg_ParseTuple(args, "OO", &py_catalogues, &py_array)) {
    return NULL;
  }

  Catalogues* const cats=
    (Catalogues*) PyCapsule_GetPointer(py_catalogues, "_Catalogues");
  py_assert_ptr(cats);

  Py_buffer view;

  //if(PyObject_GetBuffer(py_array, &view,
  //PyBUF_ANY_CONTIGUOUS | PyBUF_FORMAT) == -1) {
  if(PyObject_GetBuffer(py_array, &view, PyBUF_FORMAT | PyBUF_FULL_RO) == -1) {
    return NULL;
  }

  //
  // Decode array information
  //
  if(view.ndim != 2) {
    PyErr_SetString(PyExc_TypeError, "Expected a 2-dimensional array");
    PyBuffer_Release(&view);
    return NULL;
  }

  if(strcmp(view.format, "d") != 0) {
    PyErr_SetString(PyExc_TypeError, "Expected an array of doubles");
    PyBuffer_Release(&view);
    return NULL;
  }

  /*
  cerr << "shape[0] " << view.shape[0] << endl;
  cerr << "shape[1] " << view.shape[1] << endl;

  cerr << "stride[0] " << view.strides[0] << endl;
  cerr << "stride[1] " << view.strides[1] << endl;

  if(view.suboffsets) {
    cerr << "suboffsets " << view.suboffsets[0] << endl;
  }
  else {
    cerr << "no suboffsets\n";
  }
  */

  if(view.suboffsets) {
    PyErr_SetString(PyExc_TypeError,
		    "Array suboffsets not handled in catalogues_append");

    PyBuffer_Release(&view);
    return NULL;
  }

  
  const int n= view.shape[0];
  const int ncol= view.shape[1];

  if(n == 0) {
    PyErr_SetString(PyExc_TypeError, "no data in nbar_obs array");
    PyBuffer_Release(&view);
    return NULL;
  }
  
  if(!(ncol == 3 || ncol == 4)) {
    PyErr_SetString(PyExc_TypeError,
		    "Expected 3 or 4 columns x,y,z,weight for a catalogue");
    PyBuffer_Release(&view);
    return NULL;
  }

  double* a= (double*) view.buf;
  
  const size_t next_row= view.strides[0]/sizeof(double);
  const size_t next_col= view.strides[1]/sizeof(double);

  // New Catalogue
  Catalogue* const cat= new Catalogue();
  cat->reserve(n);
  Particle p;
  p.radec[0]= 0; p.radec[1]= 0; p.vr= 0; p.M= 0; p.flag= 0; p.w= 1;
  
  for(int i=0; i<n; ++i) {
    p.x[0]= *a;
    p.x[1]= *(a + next_col);
    p.x[2]= *(a + 2*next_col);
    if(ncol == 4)
      p.w= *(a + 3*next_col);
      
    cat->push_back(p);
  
    a += next_row;
  }

  cats->push_back(cat);

  PyBuffer_Release(&view);
  
  Py_RETURN_NONE;
}
