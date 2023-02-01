#ifndef NEUTRINOS_H_
#define NEUTRINOS_H_

#include "GRHayL.h"

// Neutrino quantities
typedef struct neutrino_luminosities {
  double nue, anue, nux;
} neutrino_luminosities;

typedef struct neutrino_opacities {
  double nue[2], anue[2], nux[2];
} neutrino_opacities;

typedef neutrino_opacities neutrino_optical_depths;

#include "NRPyLeakage.h"

#endif // NEUTRINOS_H_
