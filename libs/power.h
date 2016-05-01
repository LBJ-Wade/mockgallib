#ifndef POWER_H
#define POWER_H 1

class PowerSpectrum {
 public:
  PowerSpectrum(const char filename[]);
  ~PowerSpectrum();
  double compute_sigma(const double R) const;
  
  int n;
  double* k;
  double* P;
};

class PowerFileError {
  
};

#endif
