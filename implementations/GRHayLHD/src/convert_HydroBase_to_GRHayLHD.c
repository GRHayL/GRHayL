/********************************
 * CONVERT ET ID TO IllinoisGRMHD
 *
 * Written in 2014 by Zachariah B. Etienne
 *
 * Sets metric & MHD variables needed
 * by IllinoisGRMHD, converting from
 * HydroBase and ADMBase.
 ********************************/

#include "GRHayLHD.h"

void convert_HydroBase_to_GRHayLHD(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_convert_HydroBase_to_GRHayLHD;
  DECLARE_CCTK_PARAMETERS;

  const double poison = 0.0/0.0;
  double dummy1, dummy2, dummy3;
  double dummy4, dummy5, dummy6;

  const int imax = cctk_lsh[0];
  const int jmax = cctk_lsh[1];
  const int kmax = cctk_lsh[2];

// We use rho and press from HydroBase directly with no need to convert
#pragma omp parallel for
  for(int k=0; k<kmax; k++) {
    for(int j=0; j<jmax; j++) {
      for(int i=0; i<imax; i++) {
        const int index=CCTK_GFINDEX3D(cctkGH,i,j,k);

        rho_b[index] = rho[index];
        pressure[index] = press[index];

        const double ETvx = vel[CCTK_GFINDEX4D(cctkGH,i,j,k,0)];
        const double ETvy = vel[CCTK_GFINDEX4D(cctkGH,i,j,k,1)];
        const double ETvz = vel[CCTK_GFINDEX4D(cctkGH,i,j,k,2)];

        // IllinoisGRMHD defines v^i = u^i/u^0.

        // Meanwhile, the ET/HydroBase formalism, called the Valencia
        // formalism, splits the 4 velocity into a purely spatial part
        // and a part that is normal to the spatial hypersurface:
        // u^a = G (n^a + U^a), (Eq. 14 of arXiv:1304.5544; G=W, U^a=v^a)
        // where n^a is the unit normal vector to the spatial hypersurface,
        // n_a = {-\alpha,0,0,0}, and U^a is the purely spatial part, which
        // is defined in HydroBase as the vel[] vector gridfunction.
        // Then u^a n_a = - \alpha u^0 = G n^a n_a = -G, and
        // of course \alpha u^0 = 1/sqrt(1+γ^ij u_i u_j) = \Gamma,
        // the standard Lorentz factor.

        // Note that n^i = - \beta^i / \alpha, so
        // u^a = \Gamma (n^a + U^a)
        // -> u^i = \Gamma ( U^i - \beta^i / \alpha )
        // which implies
        // v^i = u^i/u^0
        //     = \Gamma/u^0 ( U^i - \beta^i / \alpha ) <- \Gamma = \alpha u^0
        //     = \alpha ( U^i - \beta^i / \alpha )
        //     = \alpha U^i - \beta^i

        vx[index] = alp[index]*ETvx - betax[index];
        vy[index] = alp[index]*ETvy - betay[index];
        vz[index] = alp[index]*ETvz - betaz[index];
      }
    }
  }

  // Neat feature for debugging: Add a roundoff-error perturbation
  //    to the initial data.
  // Set random_pert variable to ~1e-14 for a random 15th digit
  //    perturbation.
  if(random_pert > 1e-30) {
    srand(random_seed); // Use srand() as rand() is thread-safe.
    for(int k=0; k<kmax; k++) {
      for(int j=0; j<jmax; j++) {
        for(int i=0; i<imax; i++) {
          const int index=CCTK_GFINDEX3D(cctkGH,i,j,k);
          const double pert = (random_pert*(double)rand() / RAND_MAX);
          const double one_plus_pert=(1.0+pert);
          rho_b[index]*=one_plus_pert;
          vx[index]*=one_plus_pert;
          vy[index]*=one_plus_pert;
          vz[index]*=one_plus_pert;
        }
      }
    }
  }

  // Finally, enforce limits on primitives & compute conservative variables.
#pragma omp parallel for
  for(int k=0; k<kmax; k++) {
    for(int j=0; j<jmax; j++) {
      for(int i=0; i<imax; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        metric_quantities ADM_metric;
        ghl_initialize_metric(
              alp[index],
              betax[index], betay[index], betaz[index],
              gxx[index], gxy[index], gxz[index],
              gyy[index], gyz[index], gzz[index],
              &ADM_metric);

        ADM_aux_quantities metric_aux;
        ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

        primitive_quantities prims;
        ghl_initialize_primitives(
              rho_b[index], pressure[index], eps[index],
              vx[index], vy[index], vz[index],
              0.0, 0.0, 0.0,
              poison, poison, poison,
              &prims);

        conservative_quantities cons;
        int speed_limited = 0;
        //This applies inequality fixes on the conservatives
        ghl_enforce_primitive_limits_and_compute_u0(
              grhayl_params, grhayl_eos, &ADM_metric,
              &metric_aux, &prims, &speed_limited);
        //This computes the conservatives from the new primitives
        ghl_compute_conservs(
              &ADM_metric, &metric_aux, &prims, &cons);

        ghl_return_primitives(
              &prims,
              &rho_b[index], &pressure[index], &eps[index],
              &vx[index], &vy[index], &vz[index],
              &dummy1, &dummy2, &dummy3,
              &dummy4, &dummy5, &dummy6);

        ghl_return_conservatives(
              &cons,
              &rho_star[index], &tau[index],
              &Stildex[index], &Stildey[index], &Stildez[index],
              &dummy1, &dummy2);
      }
    }
  }
}
