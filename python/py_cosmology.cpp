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

PyObject* py_cosmology_rhom(PyObject* self, PyObject* args)
{
  const double rho_m= cosmology_rho_m();
  return Py_BuildValue("d", rho_m);
}

PyObject* py_cosmology_compute_comoving_distance(PyObject* self, PyObject* args)
{
  double a;
  if(!PyArg_ParseTuple(args, "d", &a))
    return NULL;

  return Py_BuildValue("d", cosmology_compute_comoving_distance(a));
}

