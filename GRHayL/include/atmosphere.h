#ifndef ATMOSPHERE_H_
#define ATMOSPHERE_H_

#include "grhayl.h"

#ifdef __cplusplus
extern "C" {
#endif

void grhayl_set_prims_to_constant_atm(
      const eos_parameters *restrict eos,
      primitive_quantities *restrict prims);

void grhayl_set_prims_to_radial_falloff_atm(
      const eos_parameters *restrict eos,
      const double r, //not sure what is actually needed
      primitive_quantities *restrict prims);

#ifdef __cplusplus
}
#endif

#endif // ATMOSPHERE_H
