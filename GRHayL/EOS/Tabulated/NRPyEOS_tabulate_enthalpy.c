#include "ghl_nrpyeos_tabulated.h"

// This function tabulates the specific enthalpy, assuming geometrized units.
// We store log(h), so h itself must remain strictly positive.
void NRPyEOS_tabulate_enthalpy(ghl_eos_parameters *restrict eos) {
  const double eps0 = eos->energy_shift;
  for(int iy = 0; iy < eos->N_Ye; iy++) {
    for(int it = 0; it < eos->N_T; it++) {
      for(int ir = 0; ir < eos->N_rho; ir++) {
        const int iprs = NRPYEOS_IDX3D(eos, ir, it, iy, NRPyEOS_press_key);
        const int ieps = NRPYEOS_IDX3D(eos, ir, it, iy, NRPyEOS_eps_key);
        const int ient = NRPYEOS_IDX3D(eos, ir, it, iy, NRPyEOS_enthalpy_key);

        const double rho = exp(eos->table_logrho[ir]);
        const double P = exp(eos->table_all[iprs]);
        const double eps = exp(eos->table_all[ieps]) - eps0;
        const double h = 1.0 + eps + P / rho;

        const int ih = ir + eos->N_rho * ( it + eos->N_T * iy );
        const double logh = log(h);
        eos->table_all[ient] = logh;
        eos->table_logh[ih] = logh;
      }
    }
  }
}

