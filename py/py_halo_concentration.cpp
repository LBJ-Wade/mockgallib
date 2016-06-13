#include "sigma.h"
#include "halo_concentration.h"
#include "py_assert.h"
#include "py_halo_concentration.h"



PyObject* py_halo_concentration_init(PyObject* self, PyObject* args)
{
  // halo_concentration_init(_sigma)
  PyObject* py_sigma;
  
  if(!PyArg_ParseTuple(args, "O", &py_sigma))
     return NULL;

  Sigma* const s=
    (Sigma *) PyCapsule_GetPointer(py_sigma, "_Sigma");
  py_assert_ptr(s);

  halo_concentration_init(s);

  Py_RETURN_NONE;
}


  
