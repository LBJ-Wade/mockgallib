#ifndef TINKER_H
#define TINKER_H 1

#include <gsl/gsl_integration.h>

class MF {
 public:
  MF();
  MF(const double z);

  void set_redshift(const double a);
  double f(const double nu) const;
  static double b(const double nu);
 private:
  double redshift;
  double alpha;
};

#endif

