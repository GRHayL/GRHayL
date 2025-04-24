#ifndef GHL_RADIATION_H_
#define GHL_RADIATION_H_

#include "ghl.h"

/**
 * @ingroup Radiation
 * @enum ghl_neutrino_luminosities
 * @brief  Neutrino quantities
 *
 * @todo
 * Leo, please comment this :).
 */
typedef struct ghl_neutrino_luminosities {
  double nue, anue, nux;
} ghl_neutrino_luminosities;

/**
 * @ingroup Radiation
 * @enum ghl_neutrino_luminosities
 * @brief  Neutrino quantities
 *
 * @todo
 * Leo, please comment this :).
 */
typedef struct ghl_neutrino_opacities {
  double nue[2], anue[2], nux[2];
} ghl_neutrino_opacities;

typedef ghl_neutrino_opacities ghl_neutrino_optical_depths;

#include "nrpyleakage.h"

#endif // GHL_RADIATION_H_
