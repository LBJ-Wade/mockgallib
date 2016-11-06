#ifndef CORR_PAIR_CORRECTION
#define CORR_PAIR_CORRECTION 1

void corr_pair_correction_init(const char filename[]);
void corr_pair_correction_free();
double corr_pair_correction(const double theta);

#endif
