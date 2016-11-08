#include <cassert>
#include "msg.h"
#include "comm.h"

#include "py_callback.h"
#include "py_assert.h"

#include <iostream>
using namespace std;

static int call_function_collective(PyObject* const py_function);

PyObject* py_callback_standby(PyObject *self, PyObject *args)
{
  // _callback_ready(function)
  //
  // Args:
  //    function: callable that takes a list of floats x[]
  //
  // node != 0 will call this function
  assert(comm_this_rank() != 0);
  PyObject *py_function;

  if(!PyArg_ParseTuple(args, "O", &py_function)) {
    PyErr_SetString(PyExc_TypeError, "1 arguments required for callback");
    return NULL;
  }
    
  if(!PyCallable_Check(py_function)) {
    PyErr_SetString(PyExc_TypeError, "callback function is not callable");
    return NULL;
  }
  
  while(1) {
    int command= comm_bcast_int(0);
    //cerr << "received signal " << command << endl;
    
    if(command == 123) {
      //cerr << "received signal to call callback function\n";
      int ret= call_function_collective(py_function);
      if(ret) comm_abort();
      //cerr << "call callback done\n";
    }
    else if(command == 456) {
      //cerr << "recieved signal to relase callback function\n";
      Py_RETURN_NONE;
      //break;
    }
    else
      comm_abort();
  }

  Py_RETURN_NONE;
}

PyObject* py_callback_sync(PyObject *self, PyObject *args)
{
  // _callback_exec(x)
  // called by node0 to execute callback functions in other MPI nodes
  if(comm_this_rank() != 0)
    Py_RETURN_NONE;

  assert(comm_this_rank() != 0);
  
  PyObject *py_x;
  if(!PyArg_ParseTuple(args, "O", &py_x)) {
    return NULL;
  }

  msg_printf(msg_debug, "send command 123\n");
  comm_bcast_int(123);
  
  if(!PyList_Check(py_x)) {
    PyErr_SetString(PyExc_TypeError, "x in callback_sync is not a list\n");
    return NULL;
  }
      
    
  const int n= PyList_Size(py_x);
  py_assert_ptr(n > 0);
  double* const x= (double*) malloc(sizeof(double)*n); assert(x);

  for(Py_ssize_t i=0; i<n; ++i) {
    PyObject* py_x_i= PyList_GetItem(py_x, i);
    x[i]= PyFloat_AsDouble(py_x_i);
  }

  comm_bcast_int(n);

  comm_mpi_bcast_double(x, n);
  free(x);

  Py_RETURN_NONE;
}

PyObject* py_callback_release(PyObject *self, PyObject *args)
{
  assert(comm_this_rank() == 0);
  msg_printf(msg_debug, "_callback_release 456\n");

  comm_bcast_int(456);
  Py_RETURN_NONE;
}

int call_function_collective(PyObject* const py_function)
{
  assert(comm_this_rank() != 0);

  int n= comm_bcast_int(0);
  double* const x= (double*) malloc(sizeof(double)*n);

  comm_mpi_bcast_double(x, n);

  PyObject* py_x= PyList_New(n);

  for(int i=0; i<n; ++i)
    PyList_SetItem(py_x, i, PyFloat_FromDouble(x[i]));


  // Call callback function
  PyObject* const py_arg= Py_BuildValue("(O)", py_x);
  PyObject* py_result= PyEval_CallObject(py_function, py_arg);

  Py_DECREF(py_arg);
  Py_DECREF(py_x);
  
  if(py_result == NULL) {
    msg_printf(msg_fatal, "Error occured in cost_function in minimise\n");
    PyErr_SetNone(PyExc_TypeError);
    return 1;
  }
  else if(!PyFloat_Check(py_result)) {
    msg_printf(msg_fatal,
	    "Error: return value of cost_function in minimise is not float\n");
    PyErr_SetNone(PyExc_TypeError);
    return 1;
  }

  double result= PyFloat_AsDouble(py_result);
  Py_DECREF(py_result);
  free(x);

  return 0;
}
