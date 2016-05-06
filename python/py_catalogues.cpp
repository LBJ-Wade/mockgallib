#include <stdexcept>
#include "py_assert.h"
#include "hod.h"
#include "py_catalogues.h"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"

PyMODINIT_FUNC
py_catalogues_module_init()
{
  import_array();

  return NULL;
}

static void py_catalogues_free(PyObject *obj);


PyObject* py_catalogues_alloc(PyObject* self, PyObject* args)
{
  Catalogues* const cats= new Catalogues;

  catalogue_init();

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

PyObject* py_catalogues_generate_galaxies(PyObject* self, PyObject* args)
{
  // _catalogues_generate_galaxies(_catalogues, _hod, _lightcone, z_min, z_max)
  PyObject *py_catalogues, *py_hod, *py_lightcones;
  double z_min, z_max;
  if(!PyArg_ParseTuple(args, "OOOdd", &py_catalogues, &py_hod, &py_lightcones,
		       &z_min, &z_max)) {
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

  for(size_t i=0; i<n; ++i) {
    catalogue_generate_mock(hod, lightcones->at(i),
			    z_min, z_max, cats->at(i));

    msg_printf(msg_verbose, "mock catalogue %lu generated with %lu galaxies\n",
	       i, cats->at(i)->size());
  }

  Py_RETURN_NONE;
}

PyObject* py_catalogues_len(PyObject* self, PyObject* args)
{
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
  // return ith _Catalogue
  
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
