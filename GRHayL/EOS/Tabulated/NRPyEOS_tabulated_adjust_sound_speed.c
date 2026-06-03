#include "ghl_nrpyeos_tabulated.h"

// This function performs the following tasks on the sound speed squared:
//   * Convert non-relativistic table values to relativistic values
//   * NAN          -> 1.0 - 1e-15
//   * INF          -> 1.0 - 1e-15
//   * Superluminal -> 1.0
//   * Negative     -> 0.0
void NRPyEOS_tabulated_adjust_sound_speed(ghl_eos_parameters *restrict eos, bool cs2_is_relativistic) {
  size_t nonfinite_cs2_count = 0;
  size_t superluminal_cs2_count = 0;
  size_t negative_cs2_count = 0;

  for(int iy = 0; iy < eos->N_Ye; iy++) {
    for(int it = 0; it < eos->N_T; it++) {
      for(int ir = 0; ir < eos->N_rho; ir++) {
        const int ics2 = NRPYEOS_IDX3D(eos, ir, it, iy, NRPyEOS_cs2_key);

        if(!cs2_is_relativistic) {
          const int ient = NRPYEOS_IDX3D(eos, ir, it, iy, NRPyEOS_enthalpy_key);
          const double h = exp(eos->table_all[ient]);
          eos->table_all[ics2] /= h;
        }

        if(eos->clean_sound_speed) {
          if(!isfinite(eos->table_all[ics2])) {
            nonfinite_cs2_count++;
            eos->table_all[ics2] = 1.0 - 1e-15;
          }
          else if(eos->table_all[ics2] > 1.0) {
            superluminal_cs2_count++;
            eos->table_all[ics2] = 1.0;
          }
          else if(eos->table_all[ics2] < 0.0) {
            negative_cs2_count++;
            eos->table_all[ics2] = 0.0;
          }
        }
      }
    }
  }

#define GHL_CS2_FIX_REPORT(property, fix)                                 \
  {                                                                       \
    const size_t count_ = property ## _cs2_count;                         \
    const double percent_ = 100.0 * count_ / (double)npoints;             \
    if(count_ == 0) {                                                     \
      ghl_info("No points in table have a %s cs2 value\n", #property);    \
    }                                                                     \
    else {                                                                \
      ghl_warn("Table has %zu points (~%.1f%%) with %s cs2; set to %s\n", \
               count_, percent_, #property, fix);                         \
    }                                                                     \
  }

  const int npoints = eos->N_rho * eos->N_Ye * eos->N_T;
  GHL_CS2_FIX_REPORT(nonfinite, "1.0 - 1e-15");
  GHL_CS2_FIX_REPORT(superluminal, "1.0");
  GHL_CS2_FIX_REPORT(negative, "0.0");

#undef GHL_CS2_FIX_REPORT
}
