#ifndef MF_CUMULATIVE_H
#define MF_CUMULATIVE_H 1

#include <gsl/gsl_interp.h>
#include "mf.h"
#include "sigma.h"

class MfCumulative {
 public:
  MfCumulative(const double a);
  ~MfCumulative();
  double n_cumulative(const double M) const;
  double M(const double nM) const;
  //double nM_max;
 private:
  const int n;
  double *M_array, *nM_array;
  gsl_interp *interp, *interp2;
  gsl_interp_accel *acc, *acc2;
};



#endif
