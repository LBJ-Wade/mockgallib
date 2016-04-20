//
// Fits nbar(z)
//

#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_multimin.h>

#include "msg.h"
#include "cosmology.h"
#include "nbar.h"
#include "nbar_fitting.h"

using namespace std;

static double nbar_multimin_f (const gsl_vector *v, void *params);
static double compute_cost_function(vector<Nbar> const * const vobs,
				    vector<Nbar> const * const vhod);

NbarFitting* nbar_fitting_alloc(PowerSpectrum const * const ps,
				Hod* const hod,
				const vector<Nbar>* const vnbar_obs,
				const double z_min, const double z_max)
{
  // Setup nbar fitting
  // z_min, z_max: fitting redshift range
  
  NbarFitting* fitting= new NbarFitting();

  NbarIntegration* ni0= nbar_integration_alloc(ps, hod);

  const int n= vnbar_obs->size(); assert(n > 0);
  fitting->vni.reserve(n);

  for(int i=0; i<n; ++i) {
    fitting->vni.push_back(new NbarIntegration(ni0));
  }
  assert(fitting->vni.size() == n);

  fitting->hod= hod;
  fitting->z_min= z_min;
  fitting->z_max= z_max;
  fitting->vobs= vnbar_obs;
  fitting->vhod= new vector<Nbar>(*vnbar_obs);

  msg_printf(msg_debug, "Allocated nbar_fitting\n");

  return fitting;
}

void nbar_fitting_free(NbarFitting* const fitting)
{
  const int n= fitting->vni.size();


  if(n > 0) {
    NbarIntegration* const ni= fitting->vni.front();
    nbar_integration_free(ni);
  }

  for(int i=1; i<n; ++i) {
    // i=0 is already freed above
    delete fitting->vni[i];
  }

  delete fitting->vhod;
  
  delete fitting;
}

void nbar_fitting_compute(NbarFitting* fitting)
{
  // Compute best fitting logMmin function
  //   logMmin = c[0] + c[1]*x + c[2]*x^2 + c[3]*x^3
  // for given c[4-9] coefficients in fitting->hod->c,
  // where x = z - hod::z0

  // This function updates fitting->hod->c[0-3]
  const int n= fitting->vni.size();
  
  assert(fitting->vhod->size() == n);
  assert(fitting->vobs->size() == n);

  
  const int nparam= 4;
  
  // Starting point 
  gsl_vector* const x = gsl_vector_alloc(nparam);

  gsl_vector_set(x, 0, fitting->hod->c[0]);
  gsl_vector_set(x, 1, fitting->hod->c[1]);
  gsl_vector_set(x, 2, fitting->hod->c[2]);
  gsl_vector_set(x, 3, fitting->hod->c[3]);

  // initial stepsize
  gsl_vector* const ss = gsl_vector_alloc(nparam);
  gsl_vector_set_all(ss, 0.01);
  
  // Initialize method and iterate
  gsl_multimin_function minex_func;
  minex_func.n = nparam;
  minex_func.f = nbar_multimin_f;
  minex_func.params = fitting;

  gsl_multimin_fminimizer *s=
    gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, nparam);
  gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

  int iter = 0; int status;
  const int max_iter= 200;

  
  do {
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);
    if(status)
      break;
    double size= gsl_multimin_fminimizer_size(s);
    status = gsl_multimin_test_size(size, 1e-3);


    if (status == GSL_SUCCESS) {
      fprintf(stderr, "nbar converged to minimum with %d steps.\n", iter);
    }

    msg_printf(msg_verbose, "nbar_fitting %3d %e | %.4f %.5f %.5f %.5f %.3f\n",
	   iter,
	   s->fval,
	   gsl_vector_get(s->x, 0),
	   gsl_vector_get(s->x, 1),
	   gsl_vector_get(s->x, 2),
	   gsl_vector_get(s->x, 3),
	   size);

  } while (status == GSL_CONTINUE && iter < max_iter);


  if(iter == max_iter) {
    msg_printf(msg_warn, "Reached maximum iteration for nbar fitting\n");
  }

  fitting->hod->c[0]= gsl_vector_get(s->x, 0);
  fitting->hod->c[1]= gsl_vector_get(s->x, 1);
  fitting->hod->c[2]= gsl_vector_get(s->x, 2);
  fitting->hod->c[3]= gsl_vector_get(s->x, 3);

  printf("final x %e %e %e %e\n",
	 gsl_vector_get(s->x, 0),
	 gsl_vector_get(s->x, 1),
	 gsl_vector_get(s->x, 2),
	 gsl_vector_get(s->x, 3));
  
  double chi2_final= nbar_multimin_f(x, fitting);
  printf("final %e\n", chi2_final);


  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free(s);
}


double nbar_multimin_f(const gsl_vector *v, void *params)
{
  // The minimiser library minimise this function
  NbarFitting* const fitting= (NbarFitting*) params;

  // The minisation algorithm gives the c to evaluate
  fitting->hod->c[0]= gsl_vector_get(v, 0);
  fitting->hod->c[1]= gsl_vector_get(v, 1);
  fitting->hod->c[2]= gsl_vector_get(v, 2);
  fitting->hod->c[3]= gsl_vector_get(v, 3);

  // compute n(z) from parameter hod->c[]
  const int n= fitting->vni.size();
  for(int i=0; i<n; ++i) {
    NbarIntegration* const ni= fitting->vni[i];
    const double z= (*fitting->vobs)[i].z;
    (*fitting->vhod)[i].nbar= nbar_compute(ni, z);
  }

  // evaluate the cost function between vobs and vhod
  return compute_cost_function(fitting->vobs, fitting->vhod);
}


double compute_cost_function(vector<Nbar> const * const vobs,
			     vector<Nbar> const * const vhod)
{
  // cost function is a measure how close vobs and vhod are.
  
  double chi2= 0.0;
  
  const int m= vobs->size();
  for(int i=0; i<m; ++i) {
    double diff= ((*vobs)[i].nbar - (*vhod)[i].nbar)/ (*vobs)[i].dnbar;
    chi2 += diff*diff;
    printf("%e %e %e %e\n", (*vobs)[i].z, (*vobs)[i].nbar, (*vhod)[i].nbar, diff);
  }

  return chi2;
}

