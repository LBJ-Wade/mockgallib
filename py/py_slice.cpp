#include "sky.h"
#include "remap.h"
#include "slice.h"
#include "py_assert.h"
#include "py_slice.h"

static void py_slice_free(PyObject *obj);

PyObject* py_slice_alloc(PyObject* self, PyObject* args)
{
  // _slice_alloc(_remap, _sky)
  PyObject *py_remap, *py_sky;
  if(!PyArg_ParseTuple(args, "OO", &py_remap, &py_sky))
    return NULL;

  Remap* const remap= (Remap*) PyCapsule_GetPointer(py_remap, "_Remap");
  py_assert_ptr(remap);

  Sky* const sky= (Sky*) PyCapsule_GetPointer(py_sky, "_Sky");
  py_assert_ptr(sky);


  Slice* const slice= new Slice(remap->boxsize, sky->width, sky->centre);
  
  return PyCapsule_New(slice, "_Slice", py_slice_free);
}

void py_slice_free(PyObject *obj)
{
  Slice* const slice= (Slice*) PyCapsule_GetPointer(obj, "_Slice");
  py_assert_void(slice);

  delete slice;
}

PyObject* py_slice_len(PyObject* self, PyObject* args)
{
  PyObject* py_slice;
  
  if(!PyArg_ParseTuple(args, "O", &py_slice))
     return NULL;

  Slice const * const slice= (Slice*) PyCapsule_GetPointer(py_slice, "_Slice");
  py_assert_ptr(slice);

  return Py_BuildValue("i", slice->n);
}

PyObject* py_slice_boxsize(PyObject* self, PyObject* args)
{
  PyObject* py_slice;
  
  if(!PyArg_ParseTuple(args, "O", &py_slice))
     return NULL;

  Slice const * const slice= (Slice*) PyCapsule_GetPointer(py_slice, "_Slice");
  py_assert_ptr(slice);

  return Py_BuildValue("(ddd)",
		       (double)slice->boxsize[0],
		       (double)slice->boxsize[1],
		       (double)slice->boxsize[2]);
}



