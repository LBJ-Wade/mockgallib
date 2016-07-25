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

  halo_concentration_init();

  Py_RETURN_NONE;
}


  
