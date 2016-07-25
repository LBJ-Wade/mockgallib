#ifndef RANDOM_H
#define RANDOM_H 1

void rand_init();
void rand_free();
double rand_uniform();
double rand_poisson(const double lmbda);

#endif
