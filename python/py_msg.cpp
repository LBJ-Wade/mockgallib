//
// wrapping power.cpp
//

#include "Python.h"
#include "msg.h"

PyObject* py_msg_set_loglevel(PyObject *self, PyObject* args)
{
  int level;
  if(!PyArg_ParseTuple(args, "i", &level)) {
    return NULL;
  }
  msg_set_loglevel((LogLevel) level);
  Py_RETURN_NONE;
}
