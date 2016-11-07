#include "comm.h"
#include "py_comm.h"

PyObject* py_comm_init(PyObject *self, PyObject* args)
{
  comm_init(0, 0);
  Py_RETURN_NONE;
}


PyObject* py_comm_finalise(PyObject *self, PyObject* args)
{
  comm_finalise();
  Py_RETURN_NONE;
}


PyObject* py_comm_this_rank(PyObject *self, PyObject* args)
{
  return Py_BuildValue("i", comm_this_rank());
}


PyObject* py_comm_n_nodes(PyObject *self, PyObject* args)
{
  return Py_BuildValue("i", comm_n_nodes());
}

PyObject* py_comm_bcast_str(PyObject *self, PyObject* args)
{
  char* str;
  
  if(!PyArg_ParseTuple(args, "s", &str)) {
    return NULL;
  }

  str= comm_bcast_char(str);

  PyObject* py_str= Py_BuildValue("s", str);
  free(str);

  return py_str;
}

PyObject* py_comm_bcast_int(PyObject *self, PyObject* args)
{
  int n;
  
  if(!PyArg_ParseTuple(args, "i", &n)) {
    return NULL;
  }

  n= comm_bcast_int(n);

  return Py_BuildValue("i", n);
}
