#ifndef GRHAYLIB_H_
#define GRHAYLIB_H_

#include "include/ghl.h"
#include "include/atmosphere.h"
#include "include/con2prim.h"
#include "include/induction.h"
#include "include/reconstruction.h"
#include "include/flux_source.h"
#include "include/neutrinos.h"
#include "include/nrpyeos_hybrid.h"
#include "include/nrpyeos_tabulated.h"
#include "include/nrpyleakage.h"

extern ghl_eos_parameters *ghl_eos;
extern ghl_parameters *ghl_params;
extern GRHAYL_DEVICE ghl_parameters *device_ghl_params;
extern GRHAYL_DEVICE ghl_eos_parameters *device_ghl_eos;

#endif // GRHAYLIB_H_
