//
// Information for Python
//
#include <iostream>
#include "Python.h"

#include "py_msg.h"
#include "py_power.h"

using namespace std;


static PyMethodDef methods[] = {
  {"set_loglevel", py_msg_set_loglevel, METH_VARARGS,
   "set loglevel: 0=debug, 1=verbose, ..."},
   
  {"_power_alloc", py_power_alloc, METH_VARARGS,
   "read power spectrum from file"},
  {"_power_sigma", py_power_sigma, METH_VARARGS,
   "compute sigma(R); density fluctuation smoothed with radius R"},
  {"_power_n", py_power_n, METH_VARARGS,
   "return number of power spectrum rows"},
  {"_power_k", py_power_k, METH_VARARGS,
   "return wavenumber array k"},
  {"_power_P", py_power_P, METH_VARARGS,
   "return power spectrum array P"},
  {"_power_ki", py_power_ki, METH_VARARGS,
   "return k[i]"},
  {"_power_Pi", py_power_Pi, METH_VARARGS,
   "return P[i]"},
  {NULL, NULL, 0, NULL}
};


static struct PyModuleDef module = {
  PyModuleDef_HEAD_INIT,
  "_mockgallib", // name of this module
  "A package for mock galaxy generation", // Doc String
  -1,
  methods
};

PyMODINIT_FUNC
PyInit__mockgallib(void) {
  py_power_module_init();  

  return PyModule_Create(&module);
}
