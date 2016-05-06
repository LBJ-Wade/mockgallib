#ifndef COLA_LIGHTCONE
#define COLA_LIGHTCONE

#include "snapshot.h"
#include "sky.h"
#include "remap.h"
#include "slice.h"
#include "lightcone.h"


void cola_lightcones_create(Snapshots const * const snapshots,
			    Sky const * const sky,
			    Remap const * const remap,
			    Slice const * const slice,
			    LightCones* const lightcones);


#endif
