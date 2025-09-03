#ifndef GHL_RADIATION_H_
#define GHL_RADIATION_H_

#include "ghl.h"

#ifdef __cplusplus
extern "C" {
#endif

// Neutrino quantities
typedef struct ghl_neutrino_luminosities {
  double nue, anue, nux;
} ghl_neutrino_luminosities;

typedef struct ghl_neutrino_opacities {
  double nue[2], anue[2], nux[2];
} ghl_neutrino_opacities;

typedef ghl_neutrino_opacities ghl_neutrino_optical_depths;

#ifdef __cplusplus
}
#endif

#include "ghl_nrpyleakage.h"

#endif // GHL_RADIATION_H_
