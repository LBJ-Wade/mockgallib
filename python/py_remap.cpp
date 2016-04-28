#include "py_remap.h"
#include "py_assert.h"

static void py_remap_free(PyObject *obj);

PyObject* py_remap_alloc(PyObject* self, PyObject* args)
{
  // _sigma_alloc(u, boxsize)

  PyObject* py_list;
  double boxsize;
  
  if(!PyArg_ParseTuple(args, "Od", &py_list, &boxsize))
    return NULL;

  if(!PyList_Check(py_list)) 
    PyErr_SetString(PyExc_ValueError,
		    "Expected a list for u in _remap_alloc");    
  Py_ssize_t n= PyList_Size(py_list);
  if(n != 9)
    PyErr_SetString(PyExc_ValueError,
		    "Expected 9 integers for _remap_alloc");


  int u[9];
  for(int i=0; i<9; ++i) {
    PyObject* py_val= PyList_GetItem(py_list, i);
    if(!PyLong_Check(py_val))
      PyErr_SetString(PyExc_ValueError,
		      "Expected integers in u in _remap_alloc");
    u[i] = (int) PyLong_AsLong(py_val);

    if(PyErr_Occurred())
      PyErr_SetString(PyExc_ValueError,
		      "Error occured reading u in _remap_alloc");
  }

  Remapping* const remap= (Remapping*) remap_alloc(u, boxsize);

  return PyCapsule_New(remap, "_Remapping", py_remap_free);
}

void py_remap_free(PyObject *obj)
{
  Remapping* const remap=
    (Remapping*) PyCapsule_GetPointer(obj, "_Remapping");
  py_assert(remap);

  remap_free(remap);
}

PyObject* py_remap_boxsize(PyObject* self, PyObject* args)
{
  // _remap_boxsize(_remap)
  
  PyObject* py_remap;
  
  if(!PyArg_ParseTuple(args, "O", &py_remap))
     return NULL;

  Remapping* remap= (Remapping*) PyCapsule_GetPointer(py_remap, "_Remapping");
  py_assert(remap);

  return Py_BuildValue("(ddd)",
		      remap->boxsize[0], remap->boxsize[1], remap->boxsize[2]);
}

