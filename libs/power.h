#ifndef POWER_H
#define POWER_H 1

struct PowerSpectrum {
  int n;
  double* k;
  double* P;
};

class PowerFileError {
  
};

PowerSpectrum* power_alloc(const char filename[]);
void power_free(PowerSpectrum* ps);

double power_sigma(PowerSpectrum const * const ps, const double R);

#endif
