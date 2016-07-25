#include "rand.h"

#include "py_assert.h"
#include "py_rand.h"


PyObject* py_rand_init(PyObject* self, PyObject* args)
{
  rand_init();

  Py_RETURN_NONE;
}


PyObject* py_rand_free(PyObject* self, PyObject* args)
{
  rand_free();

  Py_RETURN_NONE;
}
