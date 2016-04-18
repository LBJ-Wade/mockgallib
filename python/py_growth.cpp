//
// wrapping growth.cpp
//
#include "py_growth.h"
#include "growth.h"

PyObject* py_growth_D(PyObject* self, PyObject* args)
{
  // growth_D(a)
  
  double a;
  if(!PyArg_ParseTuple(args, "d", &a))
    return NULL;

  double D= growth_D(a);
  
  return Py_BuildValue("d", D);
}



