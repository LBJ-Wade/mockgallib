#ifndef HOD_H
#define HOD_H 1

#include <cstdio>
#include <cmath>

struct Hod {
 public:
  static const int n=12;
  Hod(const double z0_ = 0.5);

  double compute_logMmin(const double z) const {
    double x= z - z0;
    return c[0] + (c[1] + (c[2] + c[3]*x)*x)*x;
  }
  double compute_sigma(const double z) const {
    double x= z - z0;
    return c[4] + c[5]*pow(x, c[11]);
  }
  double compute_M1(const double z) const {
    double x= z - z0;
    return pow(10.0, compute_logMmin(z) + c[6] + c[7]*pow(x, c[10]));
  }
  double compute_alpha(const double z) const {
    double x= z - z0;
    return c[8] + c[9]*x;
  }
  
  void compute_param_z(const double z) {
    logMmin= compute_logMmin(z);
    sigma=   compute_sigma(z);
    M0=      pow(10.0, logMmin);
    M1=      compute_M1(z);
    alpha=   compute_alpha(z);
  }
    
  double ncen(const double M) const {
    return 0.5*(1.0 + erf((log10(M) - logMmin)/sigma));
  }
  
  double nsat(const double M) const {
    return M <= M0 ? 0.0 : pow((M - M0)/M1, alpha);
  }
  double c[n];
  const double z0;
 private:  
  double logMmin, sigma, M0, M1, alpha;
};

  

#endif
