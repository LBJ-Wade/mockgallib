//
// Mass funcation and halo bias
//

// Tinker et al (2010) ApJ 724, 878


#include <cstdio>
#include <cmath>
#include <cassert>

#include "mf.h"
#include "const.h"     // => delta_c

using namespace std;

static double integrand_mf_normalisation(double nu, void* params);

MF::MF()
{
  alpha= 0.0;
}

MF::MF(const double z)
{
  const double a= 1.0/(1.0 + z);
  set_redshift(a);
}

void MF::set_redshift(const double a)
{
  // Integral is necessary to normalise mass function (mf->alpha)
  const double z= 1.0/a - 1.0;
  if(alpha > 0.0 && redshift == z)
    return;
  
  alpha= 1.0;
  redshift= z;

  gsl_integration_cquad_workspace* w= 
    gsl_integration_cquad_workspace_alloc(100);
      
  gsl_function F;
  F.function= &integrand_mf_normalisation;
  F.params= (void*) this;
       
  double result;
  gsl_integration_cquad(&F, 1.0e-8, 10.0, 1.0e-5, 1.0e-5, w, &result, 0, 0);

  gsl_integration_cquad_workspace_free(w);

  alpha = 1.0/result;
}

double MF::f(const double nu) const
{
#ifdef DEBUG
  assert(alpha > 0.0);
#endif
  
  const double z= redshift;
  
  // Table 4 and Equations (8-12) for Delta=200
  // nu= delta_c/sigma
  // dn/dM = f(nu) rho_bar/M dln sigma^-1/dM

  //const double alpha=0.368;
  const double beta=0.589*pow(1.0 + z, 0.20);     // (9)
  const double phi= -0.729*pow(1.0 + z, -0.08);   // (10)
  const double eta= -0.243*pow(1.0 + z, 0.27);    // (11)
  const double gamma= 0.864*pow(1.0 + z, -0.01);  // (12)

  // normalization constant needs to be determined from
  // S b(nu)f(nu)dnu = 1    ... (7)
  return alpha*(1.0 + pow(beta*nu, -2.0*phi))*
         pow(nu, 2*eta)*exp(-gamma*nu*nu/2.0);    // (8)
}

double MF::b(const double nu)
{
  // assert(sigma > 0.0);
  // halo bias b(M) from Tinker et al  (2010) ApJ 724, 878

  //double nu= delta_c/sigma;

  double y= log10(200.0);
  double fy= 4.0/y; double fy2= fy*fy;
  double expfy= exp(-fy2*fy2);

  // Table 2
  double A= 1.0 + 0.24*y*expfy;                            
  double a= 0.44*y - 0.88;
  double B= 0.183;
  double b= 1.5;
  double C= 0.019 + 0.107*y + 0.19*expfy;
  double c= 2.4;
  
  return 1.0 - A*pow(nu, a)/(pow(nu, a) + pow(c::delta_c, a))
           + B*pow(nu, b) + C*pow(nu, c); // (6)
}


static double integrand_mf_normalisation(double nu, void* params)
{
  MF const * const mf= (MF const *) params;
  return mf->f(nu)*MF::b(nu);
}
