#ifndef GHL_RADIATION_H_
#define GHL_RADIATION_H_

#include "ghl.h"
#include "ghl_m1.h"
#include "ghl_neutrino_rate_provider.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @defgroup Neutrinos Neutrinos
 */

/**
 * @ingroup Neutrinos
 * @brief  Neutrino luminosities
 *
 * Luminosities for electron neutrinos, electron antineutrinos, and
 * heavy-lepton neutrinos.
 */
typedef struct ghl_neutrino_luminosities {
  double nue, anue, nux;
} ghl_neutrino_luminosities;

/**
 * @ingroup Neutrinos
 * @brief  Neutrino opacities
 *
 * Opacities for electron neutrinos, electron antineutrinos, and
 * heavy-lepton neutrinos. Each species array stores the number opacity at
 * index 0 and the energy opacity at index 1.
 */
typedef struct ghl_neutrino_opacities {
  double nue[2], anue[2], nux[2];
} ghl_neutrino_opacities;

/**
 * @ingroup Neutrinos
 * @brief Neutrino optical depths
 *
 * Uses the same species and moment indexing as ghl_neutrino_opacities,
 * but stores dimensionless optical depths rather than opacities.
 */
typedef ghl_neutrino_opacities ghl_neutrino_optical_depths;

#ifdef __cplusplus
}
#endif

#include "ghl_nrpyleakage.h"

#endif // GHL_RADIATION_H_
