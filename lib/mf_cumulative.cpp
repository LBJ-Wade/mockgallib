//
// n(>M) = \int dlnM dn/dlnM
//

#include <cassert>

#include "const.h"
#include "cosmology.h"
#include "growth.h"
#include "mf.h"
#include "sigma.h"
#include "mf_cumulative.h"


static double rho_m;
static double integrand_n_cumulative(double nu, void* params);

struct MfcParams {
  MF* mf;
  Sigma* s;
  double D;
};

MfCumulative::MfCumulative(Sigma* const s, const double a) :
  n(1001)
{
  M_array= (double*) malloc(sizeof(double)*n*2); assert(M_array);
  nM_array= M_array + n;
  
  rho_m= cosmology_rho_m();

  const double z= 1.0/a - 1.0;
  MF* const mf= mf_alloc();
  mf_set_redshift(mf, a);

  const double D= growth_D(a);
  MfcParams params;
  params.mf= mf;
  params.s=  s;
  params.D= D;

  
  gsl_integration_cquad_workspace* const w= 
    gsl_integration_cquad_workspace_alloc(100);
      
  gsl_function F;
  F.function= &integrand_n_cumulative;
  F.params= (void*) &params;

  const double logMmin= log(s->M_min);
  const double logMmax= log(s->M_max);
  const double nu_min= c::delta_c*s->sinv_min;
  const double nu_max= c::delta_c*s->sinv_max;

  for(int i=0; i<n; ++i) {
    double MM= exp(logMmin + (n-i-1)*(logMmax - logMmin)/(n-1));
    M_array[i]= MM;

    double nu= c::delta_c/D*s->sigma0_inv(MM); assert(nu >= nu_min);
    double result;
  
    gsl_integration_cquad(&F, 1.0e-8, nu_max, 1.0e-5, 1.0e-5, w, &result, 0, 0);
    nM_array[i]= result;
  }

  nM_max= nM_array[n-1];
  
  interp= gsl_interp_alloc(gsl_interp_cspline, n);
  acc= gsl_interp_accel_alloc(); 

  //gsl_interp_init(interp, M_array, nM_array, n);
  gsl_interp_init(interp, nM_array, M_array, n);

  gsl_integration_cquad_workspace_free(w);
  mf_free(mf);
}

MfCumulative::~MfCumulative()
{
  gsl_interp_free(interp);
  gsl_interp_accel_free(acc);

  free(M_array);
}

double MfCumulative::M(const double nM)
{
  return gsl_interp_eval(interp, nM_array, M_array, nM, acc);
}


double integrand_n_cumulative(double nu, void* params)
{
  MfcParams* const p= (MfcParams*) params;
  
  double sigma0= c::delta_c/(p->D*nu);
  double M= p->s->M(sigma0);

  return mf_f(p->mf, nu)*rho_m/M;
}
