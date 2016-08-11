//
// Call GSL minimizer from Python
//

#include "msg.h"
#include "py_minimise.h"
#include <gsl/gsl_multimin.h>

static PyObject* build_py_x(gsl_vector const * const x)
{
  const int n= x->size;
  PyObject* py_x= PyTuple_New(n);

  for(int i=0; i<n; ++i)
    PyTuple_SetItem(py_x, i, PyFloat_FromDouble(gsl_vector_get(x, i)));

  return py_x;
}

static int call_callback(PyObject* py_callback, gsl_vector const * const x)
{
  if(py_callback == Py_None)
    return true;
  
  PyObject* const py_x= build_py_x(x);
  PyObject* const py_arg= Py_BuildValue("(O)", py_x);
  
  PyObject* py_result= PyEval_CallObject(py_callback, py_arg);
  Py_DECREF(py_arg);
  Py_DECREF(py_x);

  if(py_result == NULL) {
    msg_printf(msg_error, "Error: calling callback function from minimise()\n");
    return false;
  }
  else {
    Py_DECREF(py_result);
  }

  return true;
}

  
static double f(const gsl_vector *x, void *params)
{
  const int n= x->size;

  PyObject* py_x= PyTuple_New(n);
  for(int i=0; i<n; ++i) {
    PyTuple_SetItem(py_x, i, PyFloat_FromDouble(gsl_vector_get(x, i)));
  }

  PyObject *py_arg= Py_BuildValue("(O)", py_x);

  PyObject* py_function= (PyObject*) params;
  PyObject* py_result= PyEval_CallObject(py_function, py_arg);
  

  Py_DECREF(py_arg);
  Py_DECREF(py_x);
  
  if(py_result == NULL || !PyFloat_Check(py_result)) {
    msg_printf(msg_fatal, "Error: return value of cost_function in minimise is not float\n");
    PyErr_SetNone(PyExc_TypeError);
    return NULL;
  }

  double result= PyFloat_AsDouble(py_result);
  Py_DECREF(py_result);


  return result;
}

PyObject* py_minimise(PyObject *self, PyObject *args)
{
  // _minimize(cost_function, callback_function, x0, ss)
  PyObject *py_cost_function, *py_callback;
  PyObject *py_x0, *py_ss;

  if(!PyArg_ParseTuple(args, "OOOO", &py_cost_function, &py_callback,
		      &py_x0, &py_ss)) {
    PyErr_SetString(PyExc_TypeError, "4 arguments required for _minimise");
    return NULL;
  }

  if (!PyCallable_Check(py_cost_function)) {
    PyErr_SetString(PyExc_TypeError, "cost_function is not callable");
    return NULL;
  }

  if (!(py_callback == Py_None || PyCallable_Check(py_callback))) {
    PyErr_SetString(PyExc_TypeError, "callback is not callable");
    return NULL;
  }

  const Py_ssize_t nparam= PyList_Size(py_x0);

  //
  // Set starting point x0
  //
  gsl_vector* const x0 = gsl_vector_alloc(nparam);
  for(Py_ssize_t i=0; i<nparam; ++i) {
    PyObject* py_x0_i= PyList_GetItem(py_x0, i);
    double x0_i= PyFloat_AsDouble(py_x0_i);
    gsl_vector_set(x0, i, x0_i);
  }

  //
  // Set stepsize ss
  //
  const Py_ssize_t nss= PyList_Size(py_x0);
  if(nss != nparam) {
      PyErr_SetString(PyExc_TypeError,
	       "The number of steps sizes is different from the number of x0");
      return NULL;
  }
  gsl_vector* const ss= gsl_vector_alloc(nparam);

  for(Py_ssize_t i=0; i<nparam; ++i) {
    PyObject* py_ss_i= PyList_GetItem(py_ss, i);
    double ss_i= PyFloat_AsDouble(py_ss_i);
    gsl_vector_set(ss, i, ss_i);
  }

  //
  // Setup minimiser
  //
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s= gsl_multimin_fminimizer_alloc(T, nparam);

  gsl_multimin_function minex_func;
  minex_func.n = nparam;
  minex_func.f = f;
  minex_func.params = py_cost_function;
  gsl_multimin_fminimizer_set(s, &minex_func, x0, ss);
  
  int iter = 0; int status;
  const int max_iter= 200;
  
  do {
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);

    if(!call_callback(py_callback, s->x))
      break;
    
    if(status) break;

    double size= gsl_multimin_fminimizer_size(s);
    status = gsl_multimin_test_size(size, 1.0e-3);
  } while (status == GSL_CONTINUE && iter < max_iter);

  if(iter == max_iter)
    msg_printf(msg_warn,
	       "minimse() reached maximum number of iteration: %d", max_iter);
  

  //
  // Build return value
  //
  PyObject* py_result= PyTuple_New(nparam);
  for(int i=0; i<nparam; ++i)
    PyTuple_SetItem(py_result, i, Py_BuildValue("d", gsl_vector_get(s->x, i)));

  //
  // Cleanup
  //
  gsl_vector_free(x0);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free(s);

  return py_result;
}
