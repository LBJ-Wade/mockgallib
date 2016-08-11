#ifndef HOD_H
#define HOD_H 1

#include <cstdio>
#include <cmath>

struct Hod {
 public:
  static const int n=10;
  Hod();
  
  void compute_param_z(const double z) {
    double x= z - z0;
    logMmin= c[0] + (c[1] + (c[2] + c[3]*x)*x)*x;
    sigma=   c[4] + c[5]*x;
    M0=      pow(10.0, logMmin);
    M1=      pow(10.0, logMmin + c[6] + c[7]*x);
    alpha=   c[8] + c[9]*x;
  }
    
  double ncen(const double M) const {
    return 0.5*(1.0 + erf((log10(M) - logMmin)/sigma));
  }
  
  double nsat(const double M) const {
    return M <= M0 ? 0.0 : pow((M - M0)/M1, alpha);
  }
  double c[n];
 private:  
  static const double z0;
  double logMmin, sigma, M0, M1, alpha;
};

  

#endif
