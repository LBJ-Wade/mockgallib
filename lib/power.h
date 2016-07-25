#ifndef POWER_H
#define POWER_H 1

void power_init(const char filename[]);
void power_free();
double power_compute_sigma(const double R);

int power_n();
double * power_P();
double * power_k();

#endif
