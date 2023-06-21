/********************************
 * CONVERT ET ID TO IllinoisGRMHD
 *
 * Written in 2014 by Zachariah B. Etienne
 *
 * Sets metric & MHD variables needed
 * by IllinoisGRMHD, converting from
 * HydroBase and ADMBase.
 ********************************/

#include "GRHayLMHD.h"

void convert_HydroBase_to_GRHayLMHD(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_convert_HydroBase_to_GRHayLMHD;
  DECLARE_CCTK_PARAMETERS;

  const int imax = cctk_lsh[0];
  const int jmax = cctk_lsh[1];
  const int kmax = cctk_lsh[2];

// We use rho and press from HydroBase directly with no need to convert
#pragma omp parallel for
  for(int k=0; k<cctk_lsh[2]; k++) {
    for(int j=0; j<cctk_lsh[1]; j++) {
      for(int i=0; i<cctk_lsh[0]; i++) {
        const int index=CCTK_GFINDEX3D(cctkGH,i,j,k);
        const int ind0=CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0);
        const int ind1=CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1);
        const int ind2=CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2);

        rho_b[index] = rho[index];
        pressure[index] = press[index];

        Ax[index] = Avec[ind0];
        Ay[index] = Avec[ind1];
        Az[index] = Avec[ind2];
        phitilde[index] = Aphi[index];

        const double ETvx = vel[ind0];
        const double ETvy = vel[ind1];
        const double ETvz = vel[ind2];

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

          phitilde[index]*=one_plus_pert;
          Ax[index]*=one_plus_pert;
          Ay[index]*=one_plus_pert;
          Az[index]*=one_plus_pert;
        }
      }
    }
  }
}
