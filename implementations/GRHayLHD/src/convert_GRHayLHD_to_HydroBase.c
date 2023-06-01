#include "GRHayLHD.h"

void convert_GRHayLHD_to_HydroBase(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_convert_GRHayLHD_to_HydroBase;
  DECLARE_CCTK_PARAMETERS;

  // Generally, we only need the HydroBase variables for diagnostic purposes, so we run the below loop only at iterations in which diagnostics are run.
  if(Convert_to_HydroBase_every==0 || cctk_iteration%Convert_to_HydroBase_every!=0) return;

#pragma omp parallel for
  for(int k=0; k<cctk_lsh[2]; k++) {
    for(int j=0; j<cctk_lsh[1]; j++) {
      for(int i=0; i<cctk_lsh[0]; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
        const int index4D0 = CCTK_GFINDEX4D(cctkGH,i,j,k,0);
        const int index4D1 = CCTK_GFINDEX4D(cctkGH,i,j,k,1);
        const int index4D2 = CCTK_GFINDEX4D(cctkGH,i,j,k,2);

        /* Note that we currently do not set Abar, Y_e, temperature, entropy, Avec[3], Aphi, Avec_stag[3], Aphi_stag */
        rho[index]   = rho_b[index];
        press[index] = pressure[index];

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
        const double lapseL=alp[index];
        const double lapseL_inv=1.0/lapseL;
        const double utU[3] = {vx[index] + betax[index],
                               vy[index] + betay[index],
                               vz[index] + betaz[index]};

        vel[index4D0] = utU[0]*lapseL_inv;
        vel[index4D1] = utU[1]*lapseL_inv;
        vel[index4D2] = utU[2]*lapseL_inv;

        // \alpha u^0 = 1/sqrt(1+γ^ij u_i u_j) = \Gamma = w_lorentz
        // First compute u^0:
        // Derivation of first equation:
        // \gamma_{ij} (v^i + \beta^i)(v^j + \beta^j)/(\alpha)^2
        //   = \gamma_{ij} 1/(u^0)^2 ( \gamma^{ik} u_k \gamma^{jl} u_l /(\alpha)^2 <- Using Eq. 53 of arXiv:astro-ph/0503420
        //   = 1/(u^0 \alpha)^2 u_j u_l \gamma^{jl}  <- Since \gamma_{ij} \gamma^{ik} = \delta^k_j
        //   = 1/(u^0 \alpha)^2 ( (u^0 \alpha)^2 - 1 ) <- Using Eq. 56 of arXiv:astro-ph/0503420
        //   = 1 - 1/(u^0 \alpha)^2 <= 1
        const double gxxL = gxx[index];
        const double gxyL = gxy[index];
        const double gxzL = gxz[index];
        const double gyyL = gyy[index];
        const double gyzL = gyz[index];
        const double gzzL = gzz[index];

        const double one_minus_one_over_alpha_u0_squared = (gxxL* SQR(utU[0]) +
                                                            2.0*gxyL*(utU[0])*(utU[1]) +
                                                            2.0*gxzL*(utU[0])*(utU[2]) +
                                                            gyyL* SQR(utU[1]) +
                                                            2.0*gyzL*(utU[1])*(utU[2]) +
                                                            gzzL* SQR(utU[2]) )*SQR(lapseL_inv);
        /*** Check for superluminal velocity ***/
        //FIXME: Instead of >1.0, should be one_minus_one_over_alpha_u0_squared > ONE_MINUS_ONE_OVER_GAMMA_SPEED_LIMIT_SQUARED, for consistency with conserv_to_prims routines

        if(one_minus_one_over_alpha_u0_squared > 1.0) {
          CCTK_VINFO("convert_from_GRHayLHD_to_HydroBase WARNING: Found superluminal velocity. This should have been caught by GRHayLHD.");
        }

        // A = 1.0-one_minus_one_over_alpha_u0_squared = 1-(1-1/(al u0)^2) = 1/(al u0)^2
        // 1/sqrt(A) = al u0
        const double alpha_u0 = 1.0/sqrt(1.0-one_minus_one_over_alpha_u0_squared);
        if(isnan(alpha_u0*lapseL_inv)) printf("BAD FOUND NAN ALPHAU0 CALC: %.15e %.15e %.15e\n",alpha_u0,lapseL_inv,one_minus_one_over_alpha_u0_squared);

        w_lorentz[index] = alpha_u0;

        Bvec[index4D0] = 0.0;
        Bvec[index4D1] = 0.0;
        Bvec[index4D2] = 0.0;
      }
    }
  }
}
