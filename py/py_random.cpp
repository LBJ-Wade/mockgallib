#include <gsl/gsl_rng.h>
#include <time.h> 

#include "py_assert.h"
#include "py_random.h"

static void py_random_free(PyObject *obj);
  
PyObject* py_random_alloc(PyObject* self, PyObject* args)
{
  const unsigned int seed= (unsigned int) time(NULL);
  gsl_rng* const rng= gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(rng, seed);
 
  return PyCapsule_New(rng, "_GSL_RNG", py_random_free);
}

void py_random_free(PyObject *obj)
{
  gsl_rng* const rng= (gsl_rng*) PyCapsule_GetPointer(obj, "_GSL_RNG");

  py_assert_void(rng);

  gsl_rng_free(rng);
}


