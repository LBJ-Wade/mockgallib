//
// Call GSL minimizer from Python
//

#include <iostream>
#include <cstdlib>
#include <gsl/gsl_multimin.h>

#include "msg.h"
#include "comm.h"
#include "py_minimise.h"

using namespace std;
enum MinimizeCommand {min_evaluate_f=0, min_log=1, min_end=2};

static PyObject* build_py_x_collective(gsl_vector const * const x);

static double call_cost_function_collective(PyObject* const py_cost_function,
					    const gsl_vector *x);
static int call_log_function_collective(PyObject* const py_log_function,
					gsl_vector const * const x);

static PyObject* root_minimiser(PyObject * const py_cost_function,
			 PyObject * const py_callback,
			 PyObject * const py_x0,
			 PyObject * const py_ss);
static void minimiser_loop(PyObject* const py_cost_function,
			   PyObject* const py_log_function);

double f(const gsl_vector *x, void *params);

//
// Main minimise function
//
PyObject* py_minimise(PyObject *self, PyObject *args)
{
  // _minimize(cost_function, log_function, x0, ss)
  //
  // Args:
  //    log_function: callable or None
  // Returns:
  //    Tuple of x in rank = 0
  //    None in rank > 0

  PyObject *py_cost_function, *py_log_function;
  PyObject *py_x0, *py_ss;

  if(!PyArg_ParseTuple(args, "OOOO", &py_cost_function, &py_log_function,
		       &py_x0, &py_ss)) {
    PyErr_SetString(PyExc_TypeError, "4 arguments required for _minimise");
    return NULL;
  }
    
  if (!PyCallable_Check(py_cost_function)) {
    PyErr_SetString(PyExc_TypeError, "cost_function is not callable");
    return NULL;
  }
  
  if (!(py_log_function == Py_None || PyCallable_Check(py_log_function))) {
    PyErr_SetString(PyExc_TypeError, "log_funtion is not callable");
    return NULL;
  }

  if(comm_this_rank() == 0) {
     PyObject* py_result=
       root_minimiser(py_cost_function, py_log_function, py_x0, py_ss);
     return py_result;
  }
  else
    minimiser_loop(py_cost_function, py_log_function);

  Py_RETURN_NONE;
}

PyObject* root_minimiser(PyObject * const py_cost_function,
		    PyObject * const py_log_function,
		    PyObject * const py_x0,
		    PyObject * const py_ss)
{
  // Iterate GSL minimiser in node = 0
  assert(comm_this_rank() == 0);

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
  gsl_multimin_fminimizer * const s= gsl_multimin_fminimizer_alloc(T, nparam);

  gsl_multimin_function minex_func;
  minex_func.n= nparam;
  minex_func.f= f;
  minex_func.params= py_cost_function;
  gsl_multimin_fminimizer_set(s, &minex_func, x0, ss);

  //
  // Minimisation loop
  //
  int iter = 0; int status;
  const int max_iter= 200;
  int log_status;
  
  do {
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);

    comm_bcast_int(min_log);
    log_status= call_log_function_collective(py_log_function, s->x);
      
    //if(log_status || status) break;
    if(log_status == false) {
      msg_printf(msg_fatal, "log_status error\n");
      break;
    }

    double size= gsl_multimin_fminimizer_size(s);
    msg_printf(msg_debug, "gsl minimiser size %.4f\n", size);
    status = gsl_multimin_test_size(size, 1.0e-3);
  } while (status == GSL_CONTINUE && iter < max_iter);

  if(iter == max_iter)
    msg_printf(msg_warn,
	       "minimse() reached maximum number of iteration: %d", max_iter);
  

  //
  // Build return value
  //
  PyObject* py_result= NULL;
  
  if(log_status) {
    py_result= PyTuple_New(nparam);
    for(int i=0; i<nparam; ++i)
      PyTuple_SetItem(py_result, i,
		      Py_BuildValue("d", gsl_vector_get(s->x, i)));
  }
  else {
    PyErr_SetString(PyExc_TypeError,
		    "Error: calling log function from minimise()");
  }

  //
  // Cleanup
  //
  gsl_vector_free(x0);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free(s);

  comm_bcast_int(min_end);

  return py_result;
}

void minimiser_loop(PyObject* const py_cost_function,
		    PyObject* const py_log_function)
{
  assert(comm_this_rank() != 0);
  
  while(1) {
    int command= comm_bcast_int(0);
    
    switch(command) {
    case min_evaluate_f:
      call_cost_function_collective(py_cost_function, NULL);
      break;
      
    case min_log:
      call_log_function_collective(py_log_function, NULL);
      break;
      
    case min_end:
      return;
    }
  }
}


//
// static functions
//
PyObject* build_py_x_collective(gsl_vector const * const x)
{
  // All MPI node returns a Python Tuple that corresponds to a
  // GSL vector x in node 0
  //
  // Collective call (must be called by all MPI nodes)
  // Args x: gsl_vector in node 0
  //       : not used in node > 0

  
  int n= 0;
  if(comm_this_rank() == 0)
    n= x->size;
  
  n= comm_bcast_int(n);
  PyObject* py_x= PyTuple_New(n);

  double* const xx= (double*) malloc(sizeof(double)*n); assert(xx);

  if(comm_this_rank() == 0) {
    for(int i=0; i<n; ++i)
      xx[i]= gsl_vector_get(x, i);
  }
  comm_mpi_bcast_double(xx, n);

  for(int i=0; i<n; ++i)
    PyTuple_SetItem(py_x, i, PyFloat_FromDouble(xx[i]));

  free(xx);

  return py_x;
}


double call_cost_function_collective(PyObject* const py_cost_function,
				     const gsl_vector *x)
				     
{
  // Evaluate cost function
  PyObject* const py_x= build_py_x_collective(x);
  PyObject* const py_arg= Py_BuildValue("(O)", py_x);

  PyObject* py_result= PyEval_CallObject(py_cost_function, py_arg);

  Py_DECREF(py_arg);
  Py_DECREF(py_x);
  
  if(py_result == NULL) {
    msg_printf(msg_fatal, "Error occured in cost_function in minimise\n");
    PyErr_SetNone(PyExc_TypeError);
    return 0.0;
  }
  else if(!PyFloat_Check(py_result)) {
    msg_printf(msg_fatal,
	    "Error: return value of cost_function in minimise is not float\n");
    PyErr_SetNone(PyExc_TypeError);
    return 0.0;
  }

  double result= PyFloat_AsDouble(py_result);
  Py_DECREF(py_result);

  return result;
}

 
double f(const gsl_vector *x, void *params)
{
  // Function called from GSL minimiser; the function minimised
  PyObject* py_cost_function= (PyObject*) params;

  comm_bcast_int(min_evaluate_f);
  return call_cost_function_collective(py_cost_function, x);
}


int call_log_function_collective(PyObject* const py_log_function,
				 gsl_vector const * const x)
{
  //
  // Call Python function py_log_function (collective)
  //
  if(py_log_function == Py_None)
    return true;
  
  PyObject* const py_x= build_py_x_collective(x);
  PyObject* const py_arg= Py_BuildValue("(O)", py_x);
  
  PyObject* py_result= PyEval_CallObject(py_log_function, py_arg);
  Py_DECREF(py_arg);
  Py_DECREF(py_x);

  int int_result;
  if(py_result == NULL) {
    msg_printf(msg_error, "Error: calling callback function from minimise()\n");
    int_result= 0;
  }
  else {
    Py_DECREF(py_result);
    int_result= 1;
  }

  int_result= comm_allreduce_min_int(int_result);

  if(int_result == 0)
    PyErr_SetString(PyExc_TypeError,
		    "Error: calling log function from minimise()");

  
  return int_result;
}

