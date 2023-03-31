/********************************
 * CONVERT ET ID TO IllinoisGRMHD
 *
 * Written in 2014 by Zachariah B. Etienne
 *
 * Sets metric & MHD variables needed
 * by IllinoisGRMHD, converting from
 * HydroBase and ADMBase.
 ********************************/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "IGH.h"

void convert_from_ADMBase_HydroBase_to_GRHayL_IGH(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_convert_from_ADMBase_HydroBase_to_GRHayL_IGH;
  DECLARE_CCTK_PARAMETERS;

  const double poison = 0.0/0.0;
  double dummy1, dummy2, dummy3;
  double dummy4, dummy5, dummy6;

  // Convert ADM variables (from ADMBase) to the BSSN-based variables expected by this routine.
  GRHayL_IGH_convert_ADM_to_BSSN(cctkGH,
                             gxx, gxy, gxz, gyy, gyz, gzz,
                             phi_bssn, psi_bssn,
                             gtxx, gtxy, gtxz, gtyy, gtyz, gtzz,
                             gtupxx, gtupxy, gtupxz, gtupyy, gtupyz, gtupzz);

  const int imax = cctk_lsh[0];
  const int jmax = cctk_lsh[1];
  const int kmax = cctk_lsh[2];

// We use rho and press from HydroBase directly with no need to convert
#pragma omp parallel for
  for(int k=0; k<kmax; k++)
    for(int j=0; j<jmax; j++)
      for(int i=0; i<imax; i++) {
        const int index=CCTK_GFINDEX3D(cctkGH,i,j,k);

        //TODO: this references simple gamma law; clearly needs to be extended to more EOS types
        // P = (\Gamma - 1) rho epsilon
        // -> \Gamma = P/(rho epsilon) + 1
        const double measured_gamma = ( press[index]/(rho[index] * eps[index]) + 1.0 );
        if(rho[index]>grhayl_eos->rho_atm && fabs(grhayl_eos->Gamma_th - measured_gamma)/grhayl_eos->Gamma_th > 1e-2)
          CCTK_VERROR("Expected simple gamma law with gamma_th=%.15e, but found a point with gamma law such that gamma_th=%.15e. error = %e| rb=%e rbatm=%e P=%e\n",
                      grhayl_eos->Gamma_th, measured_gamma, (grhayl_eos->Gamma_th-measured_gamma)/grhayl_eos->Gamma_th, rho[index], grhayl_eos->rho_atm, press[index] );

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

  // Neat feature for debugging: Add a roundoff-error perturbation
  //    to the initial data.
  // Set random_pert variable to ~1e-14 for a random 15th digit
  //    perturbation.
  srand(random_seed); // Use srand() as rand() is thread-safe.
  for(int k=0; k<kmax; k++)
    for(int j=0; j<jmax; j++)
      for(int i=0; i<imax; i++) {
        const int index=CCTK_GFINDEX3D(cctkGH,i,j,k);
        const double pert = (random_pert*(double)rand() / RAND_MAX);
        const double one_plus_pert=(1.0+pert);
        rho_b[index]*=one_plus_pert;
        vx[index]*=one_plus_pert;
        vy[index]*=one_plus_pert;
        vz[index]*=one_plus_pert;
  }

//TODO: comment not applicable anymore
  // FIXME: IllinoisGRMHD's Conservative-to-Primitive solver
  //   (a.k.a., C2P or con2prim) only implements single gamma-law EOS.
  // Also, not compatible with EOS driver in ET, so EOS parameters must
  //   be specified for both initial data thorns AND IllinoisGRMHD
  // TODO: Incorporate checks to ensure compatibility with ID.
  //   Alternatively, read in EOS stuff from an ET EOS driver.

  // Finally, enforce limits on primitives & compute conservative variables.
#pragma omp parallel for
  for(int k=0; k<kmax; k++)
    for(int j=0; j<jmax; j++)
      for(int i=0; i<imax; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        metric_quantities metric;
        initialize_metric(alp[index], gxx[index], gxy[index], gxz[index],
                          gyy[index], gyz[index], gzz[index], betax[index],
                          betay[index], betaz[index], &metric);

        primitive_quantities prims;
        initialize_primitives(
                          rho_b[index], pressure[index], eps[index],
                          vx[index], vy[index], vz[index],
                          0.0, 0.0, 0.0,
                          poison, poison, poison,
                          &prims);
//TODO: add support for other vars; might need to depend on whether these vars
//      have storage allocated by looking at params
//                          entropy[index], Y_e[index], temperature[index],

        conservative_quantities cons;
        stress_energy Tmunu;
        int speed_limited = 0;
        //This applies the inequality (or "Faber") fixes on the conservatives
        enforce_primitive_limits_and_compute_u0(grhayl_params, grhayl_eos, &metric, &prims, &speed_limited);
        //This computes the conservatives and stress-energy tensor from the new primitives
        compute_conservs_and_Tmunu(grhayl_params, grhayl_eos, &metric, &prims, &cons, &Tmunu);

        return_primitives(&prims,
                          &rho_b[index], &pressure[index], &eps[index],
                          &vx[index], &vy[index], &vz[index],
                          &dummy1, &dummy2, &dummy3,
                          &dummy4, &dummy5, &dummy6);
//TODO: add support for other vars; might need to depend on whether these vars
//      have storage allocated by looking at params
//                          &entropy[index], &Y_e[index], &temperature[index]);

        return_conservatives(&cons,
                          &rho_star[index], &tau[index],
                          &Stildex[index], &Stildey[index], &Stildez[index],
                          &dummy1, &dummy2);
//TODO: add support for other vars; might need to depend on whether these vars
//      have storage allocated by looking at params
//                          &Y_e[index], &entropy[index]);

        if(grhayl_params->update_Tmunu) {
          eTtt[index] = Tmunu.Ttt;
          eTtx[index] = Tmunu.Ttx;
          eTty[index] = Tmunu.Tty;
          eTtz[index] = Tmunu.Ttz;
          eTxx[index] = Tmunu.Txx;
          eTxy[index] = Tmunu.Txy;
          eTxz[index] = Tmunu.Txz;
          eTyy[index] = Tmunu.Tyy;
          eTyz[index] = Tmunu.Tyz;
          eTzz[index] = Tmunu.Tzz;
        }
  }
}
