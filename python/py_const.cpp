//
// wrapping power.cpp
//
#include "Python.h"
#include "const.h"

using namespace std;

PyObject* py_const_G(PyObject* self, PyObject* args)
{
  return Py_BuildValue("d", G);
}

PyObject* py_const_rhocrit0(PyObject* self, PyObject* args)
{
  return Py_BuildValue("d", rho_crit_0);
}

