#ifndef GRHAYLIB_H_
#define GRHAYLIB_H_

#include "include/ghl.h"
#include "include/ghl_atmosphere.h"
#include "include/ghl_con2prim.h"
#include "include/ghl_induction.h"
#include "include/ghl_reconstruction.h"
#include "include/ghl_flux_source.h"
#include "include/ghl_radiation.h"
#include "include/nrpyeos_hybrid.h"
#include "include/nrpyeos_tabulated.h"
#include "include/nrpyleakage.h"

extern ghl_eos_parameters *ghl_eos;
extern ghl_parameters *ghl_params;

#endif // GRHAYLIB_H_
