#include "py_distance.h"

PyObject* py_distance_init(PyObject* self, PyObject* args)
{
  // _distance_init(z_max)

  double z_max;
  
  if(!PyArg_ParseTuple(args, "d", &z_max))
    return NULL;
  
  distance_init(z_max);

  return Py_BuildValue("d", distance_max());
}

PyObject* py_distance_redshift(PyObject* self, PyObject* args)
{
  // _distance_redshift(d) = z(d)

  double d;
  if(!PyArg_ParseTuple(args, "d", &d))
    return NULL;
    

  return Py_BuildValue("d", distance_redshift(d));
}
