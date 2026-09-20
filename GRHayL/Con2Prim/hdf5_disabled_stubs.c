#include "ghl_con2prim.h"

#ifdef GHL_DISABLE_HDF5

#define GHL_DISABLED_TABULATED_C2P(name)                                               \
  ghl_error_codes_t name(                                                              \
        const ghl_parameters *restrict params, const ghl_eos_parameters *restrict eos, \
        const ghl_metric_quantities *restrict metric_adm,                              \
        const ghl_ADM_aux_quantities *restrict metric_aux,                             \
        const ghl_conservative_quantities *restrict cons,                              \
        ghl_primitive_quantities *restrict prims,                                      \
        ghl_con2prim_diagnostics *restrict diagnostics) {                              \
    (void)params;                                                                      \
    (void)eos;                                                                         \
    (void)metric_adm;                                                                  \
    (void)metric_aux;                                                                  \
    (void)cons;                                                                        \
    (void)prims;                                                                       \
    (void)diagnostics;                                                                 \
    return ghl_error_used_disabled_hdf5;                                               \
  }

GHL_DISABLED_TABULATED_C2P(ghl_tabulated_Noble2D)
GHL_DISABLED_TABULATED_C2P(ghl_tabulated_Palenzuela1D_energy)
GHL_DISABLED_TABULATED_C2P(ghl_tabulated_Palenzuela1D_entropy)
GHL_DISABLED_TABULATED_C2P(ghl_tabulated_Newman1D_energy)
GHL_DISABLED_TABULATED_C2P(ghl_tabulated_Newman1D_entropy)

#else
typedef int ghl_hdf5_disabled_stubs_translation_unit_is_not_empty;
#endif
