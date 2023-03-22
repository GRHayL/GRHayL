#include "cctk.h"
#include "IGH.h"

static double eos_Gamma_eff(const eos_parameters *restrict eos, const double rho_in, const double press_in);

void GRHayL_IGH_reconstruction_loop(const cGH *restrict cctkGH, const int flux_dir, const int num_vars, 
                         const int *restrict var_indices,
                         const eos_parameters *restrict eos,
                         const double **in_prims,
                         double **out_prims_r,
                         double **out_prims_l) {

  // Bounds are determined by the stencil, which requires a ghostzone of at least
  // 3, but upper index includes first ghostzone point (stencil is only 2 on upper end)
  // This limit only applies to the direction of the stencil, hence the == logic below.

  const int xdir = (flux_dir == 0);
  const int ydir = (flux_dir == 1);
  const int zdir = (flux_dir == 2);

  const int imin = cctkGH->cctk_nghostzones[0]*xdir;
  const int imax = cctkGH->cctk_lsh[0] - (cctkGH->cctk_nghostzones[0]-1)*xdir;
  const int jmin = cctkGH->cctk_nghostzones[1]*ydir;
  const int jmax = cctkGH->cctk_lsh[1] - (cctkGH->cctk_nghostzones[1]-1)*ydir;
  const int kmin = cctkGH->cctk_nghostzones[2]*zdir;
  const int kmax = cctkGH->cctk_lsh[2] - (cctkGH->cctk_nghostzones[2]-1)*zdir;

#pragma omp parallel for
  for(int k=kmin; k<kmax; k++)
    for(int j=jmin; j<jmax; j++)
      for(int i=imin; i<imax; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH, i, j, k);
        double rho[6], pressure[6], v_flux_dir[6];
        double rhor, rhol, pressr, pressl;
        double var_data[num_vars][6], vars_r[num_vars], vars_l[num_vars];

        for(int ind=0; ind<6; ind++) {
          const int stencil = CCTK_GFINDEX3D(cctkGH, i+xdir*(ind-3), j+ydir*(ind-3), k+zdir*(ind-3)); // PPM needs indices from -3 to +2
          v_flux_dir[ind] = in_prims[VX+flux_dir][stencil]; // Could be smaller; doesn't use full stencil
          rho[ind] = in_prims[RHOB][stencil];
          pressure[ind] = in_prims[PRESSURE][stencil];
          for(int var=0; var<num_vars; var++) {
            var_data[var][ind] = in_prims[var_indices[var]][stencil];
          }
        }

        // Compute Gamma
        const double Gamma = eos_Gamma_eff(eos, in_prims[RHOB][index], in_prims[PRESSURE][index]);

        simple_ppm(
          rho, pressure, var_data, num_vars,
          v_flux_dir, Gamma,
          &rhor, &rhol, &pressr, &pressl, vars_r, vars_l);

        out_prims_r[RHOB][index] = rhor;
        out_prims_l[RHOB][index] = rhol;
        out_prims_r[PRESSURE][index] = pressr;
        out_prims_l[PRESSURE][index] = pressl;
        for(int var=0; var<num_vars; var++) {
          out_prims_r[var_indices[var]][index] = vars_r[var];
          out_prims_l[var_indices[var]][index] = vars_l[var];
        }
  }
}

static double eos_Gamma_eff(const eos_parameters *restrict eos, const double rho_in, const double press_in) {
  double K, Gamma;
  eos->hybrid_get_K_and_Gamma(eos, rho_in, &K, &Gamma);
  const double P_cold = K*pow(rho_in, Gamma);
  return eos->Gamma_th + (Gamma - eos->Gamma_th)*P_cold/press_in;
}
