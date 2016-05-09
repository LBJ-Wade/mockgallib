#ifndef COLA_FILE_H
#define COLA_FILE_H 1

#include <vector>
#include "halo.h"

void cola_halo_file_open(const char filename[], float* const boxsize);
int  cola_halo_file_read_one(Halo* const h);
void cola_halo_file_close();

#endif
