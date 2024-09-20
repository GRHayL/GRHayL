#ifndef GHL_ATMOSPHERE_H_
#define GHL_ATMOSPHERE_H_

#include "ghl.h"

#ifdef __cplusplus
extern "C" {
#endif

void ghl_set_prims_to_constant_atm(
      const ghl_eos_parameters *restrict eos,
      ghl_primitive_quantities *restrict prims);

void ghl_set_prims_to_radial_falloff_atm(
      const ghl_eos_parameters *restrict eos,
      const double r, //not sure what is actually needed
      ghl_primitive_quantities *restrict prims);

#ifdef __cplusplus
}
#endif

#endif // GHL_ATMOSPHERE_H
