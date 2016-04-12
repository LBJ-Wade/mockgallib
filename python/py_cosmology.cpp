//
// Wrapping cosmology.cpp
//

#include "Python.h"
#include "cosmology.h"

PyObject* py_cosmology_set(PyObject* self, PyObject* args)
{
  // py_cosmology_set(omega_m)
  double omega_m;
  if(!PyArg_ParseTuple(args, "d", &omega_m))
    return NULL;

  cosmology_set_omega_m(omega_m);

  Py_RETURN_NONE;
}

