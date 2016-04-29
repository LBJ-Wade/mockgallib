#include "lightcone.h"
#include "remap.h"
#include "slice.h"
#include "cola_file.h"

LightCones::~LightCones()
{
  for(LightCones::iterator p= begin(); p != end(); ++p) {
    delete *p;
  }
}

void fill_lightcone(const char filename[],
		      Remapping const * const remap,
		      Slice const * const slice)
{
  cola_halo_file_open(filename);

  cola_halo_file_close();
