#ifndef PY_ASSERT_H
#define PY_ASSERT_H 1

#include "Python.h"
#include "msg.h"

// borrowed from <assert.h>
#define py_assert(e)  \
    ((void) ((e) ? ((void)0) : __py_assert (#e, __FILE__, __LINE__)))

#define __py_assert(e, file, line) \
  ((void)msg_printf(msg_fatal, "%s:%u: failed assertion `%s'\n", file, line, e), PyErr_SetString(PyExc_TypeError, "py_assersion error"))

#endif
