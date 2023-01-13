#include "cctk.h"
#include "cctk_Parameters.h"
#include "NRPyLeakageET.h"

#define velx (&vel[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely (&vel[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz (&vel[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])

#define CHECK_POINTER(pointer,name) \
  if( !pointer ) CCTK_VERROR("Failed to get pointer for gridfunction '%s'",name);

void NRPyLeakageET_compute_neutrino_opacities_and_add_source_terms_to_MHD_rhss(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(verbosity_level>1) CCTK_VINFO("Inside NRPyLeakageET_compute_opacities_and_add_source_terms_to_MHD_rhss");

  const int timelevel = 0;

  // Step 1: Get pointers to opacity and optical depth gridfunctions
  CCTK_REAL *Ye_star_rhs = (CCTK_REAL *)(CCTK_VarDataPtr(cctkGH,timelevel,GFstring_Ye_star_rhs));
  CCTK_REAL *tau_rhs     = (CCTK_REAL *)(CCTK_VarDataPtr(cctkGH,timelevel,GFstring_tau_rhs));
  CCTK_REAL *st_x_rhs    = (CCTK_REAL *)(CCTK_VarDataPtr(cctkGH,timelevel,GFstring_st_x_rhs));
  CCTK_REAL *st_y_rhs    = (CCTK_REAL *)(CCTK_VarDataPtr(cctkGH,timelevel,GFstring_st_y_rhs));
  CCTK_REAL *st_z_rhs    = (CCTK_REAL *)(CCTK_VarDataPtr(cctkGH,timelevel,GFstring_st_z_rhs));

  // Step 2: Check pointers are ok
  CHECK_POINTER(Ye_star_rhs, GFstring_Ye_star_rhs);
  CHECK_POINTER(tau_rhs    , GFstring_tau_rhs    );
  CHECK_POINTER(st_x_rhs   , GFstring_st_x_rhs   );
  CHECK_POINTER(st_y_rhs   , GFstring_st_y_rhs   );
  CHECK_POINTER(st_z_rhs   , GFstring_st_z_rhs   );

  // Step 3: Ghostzones begin and end index
  const int imin = cctk_nghostzones[0];
  const int imax = cctk_lsh[0] - cctk_nghostzones[0];
  const int jmin = cctk_nghostzones[1];
  const int jmax = cctk_lsh[1] - cctk_nghostzones[1];
  const int kmin = cctk_nghostzones[2];
  const int kmax = cctk_lsh[2] - cctk_nghostzones[2];


  // Step 3: Compute opacities and leakage source terms
  int nan_found=0,num_points=0;
  CCTK_REAL Ye_star_rhs_avg=0,tau_rhs_avg=0,st_x_rhs_avg=0,st_y_rhs_avg=0,st_z_rhs_avg=0;
#pragma omp parallel for reduction(+:nan_found,num_points,Ye_star_rhs_avg,tau_rhs_avg,st_x_rhs_avg,st_y_rhs_avg,st_z_rhs_avg)
  for(int k=kmin;k<kmax;k++) {
    for(int j=jmin;j<jmax;j++) {
      for(int i=imin;i<imax;i++) {

        // Step 3.a: Set the index of the current gridpoint
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        // Step 3.b: Check if we are within the threshold
        const CCTK_REAL rhoL = rho[index];
        if( rhoL < rho_min_threshold || rhoL > rho_max_threshold ) {
          // Step 3.b.i: Below density threshold.
          //             Set opacities to zero; don't add anything to the RHSs
          kappa_0_nue [index] = 0.0;
          kappa_1_nue [index] = 0.0;
          kappa_0_anue[index] = 0.0;
          kappa_1_anue[index] = 0.0;
          kappa_0_nux [index] = 0.0;
          kappa_1_nux [index] = 0.0;
        }
        else {
          CCTK_REAL gxxL        = gxx[index];
          CCTK_REAL gxyL        = gxy[index];
          CCTK_REAL gxzL        = gxz[index];
          CCTK_REAL gyyL        = gyy[index];
          CCTK_REAL gyzL        = gyz[index];
          CCTK_REAL gzzL        = gzz[index];
          const CCTK_REAL gdet  = fabs(gxxL * gyyL * gzzL
                                     + gxyL * gyzL * gxzL
                                     + gxzL * gxyL * gyzL
                                     - gxzL * gyyL * gxzL
                                     - gxyL * gxyL * gzzL
                                     - gxxL * gyzL * gyzL);
          const CCTK_REAL phiL  = (1.0/12.0) * log(gdet);
          const CCTK_REAL psiL  = exp(phiL);
          const CCTK_REAL psi2L = psiL *psiL;
          const CCTK_REAL psi4L = psi2L*psi2L;
          const CCTK_REAL psi6L = psi4L*psi2L;
          if( psi6L > psi6_threshold ) {
            kappa_0_nue [index] = 0.0;
            kappa_1_nue [index] = 0.0;
            kappa_0_anue[index] = 0.0;
            kappa_1_anue[index] = 0.0;
            kappa_0_nux [index] = 0.0;
            kappa_1_nux [index] = 0.0;
          }
          else {
            // Step 3.c: Read from main memory
            const CCTK_REAL alpL         = alp[index];
            const CCTK_REAL alpinvsqrdL  = 1.0/(alpL*alpL);
            const CCTK_REAL betaxL       = betax[index];
            const CCTK_REAL betayL       = betay[index];
            const CCTK_REAL betazL       = betaz[index];
            CCTK_REAL vxL                = alpL*velx[index] - betaxL;
            CCTK_REAL vyL                = alpL*vely[index] - betayL;
            CCTK_REAL vzL                = alpL*velz[index] - betazL;
            const CCTK_REAL Y_eL         = Y_e[index];
            const CCTK_REAL temperatureL = temperature[index];
            neutrino_optical_depths tauL;
            tauL.nue [0] = tau_0_nue_p [index];
            tauL.nue [1] = tau_1_nue_p [index];
            tauL.anue[0] = tau_0_anue_p[index];
            tauL.anue[1] = tau_1_anue_p[index];
            tauL.nux [0] = tau_0_nux_p [index];
            tauL.nux [1] = tau_1_nux_p [index];           

            // Step 3.d: Compute BSSN quantities; enforce det(gammabar_{ij}) = 1
            // Step 3.d.i: Compute the determinant of the physical metric

            // Step 3.d.ii: Compute BSSN quantities (gb = "gbar" = conformal metric)
            const CCTK_REAL psim4L = 1.0/(psi4L);
            CCTK_REAL gbxxL        = gxxL*psim4L;
            CCTK_REAL gbxyL        = gxyL*psim4L;
            CCTK_REAL gbxzL        = gxzL*psim4L;
            CCTK_REAL gbyyL        = gyyL*psim4L;
            CCTK_REAL gbyzL        = gyzL*psim4L;
            CCTK_REAL gbzzL        = gzzL*psim4L;
            // Step 3.d.iii: Compute the determinant of the conformal metric
            CCTK_REAL gbdet = fabs(gbxxL * gbyyL * gbzzL
                                 + gbxyL * gbyzL * gbxzL
                                 + gbxzL * gbxyL * gbyzL
                                 - gbxzL * gbyyL * gbxzL
                                 - gbxyL * gbxyL * gbzzL
                                 - gbxxL * gbyzL * gbyzL);

            // Step 3.d.iv: Enforce det(gammabar) = 1 constraint
            CCTK_REAL gbdet_Fm1o3 = fabs(1.0/cbrt(gbdet));
            gbxxL *= gbdet_Fm1o3;
            gbxxL *= gbdet_Fm1o3;
            gbxxL *= gbdet_Fm1o3;
            gbxxL *= gbdet_Fm1o3;
            gbxxL *= gbdet_Fm1o3;
            gbxxL *= gbdet_Fm1o3;
            // Step 3.d.v: Recompute physical metric
            gxxL = gbxxL*psi4L;
            gxyL = gbxyL*psi4L;
            gxzL = gbxzL*psi4L;
            gyyL = gbyyL*psi4L;
            gyzL = gbyzL*psi4L;
            gzzL = gbzzL*psi4L;

            // Step 3.e: Compute u^{mu}
            // Step 3.e.i: Compute Lorentz factor W
            const CCTK_REAL one_minus_one_over_Wsqrd = (    gxxL*(vxL + betaxL)*(vxL + betaxL) +
                                                        2.0*gxyL*(vxL + betaxL)*(vyL + betayL) +
                                                        2.0*gxzL*(vxL + betaxL)*(vzL + betazL) +
                                                            gyyL*(vyL + betayL)*(vyL + betayL) +
                                                        2.0*gyzL*(vyL + betayL)*(vzL + betazL) +
                                                            gzzL*(vzL + betazL)*(vzL + betazL) )*alpinvsqrdL;
            CCTK_REAL W = 1.0/sqrt( 1 - one_minus_one_over_Wsqrd );

            // Step 3.e.ii: Impose a speed limit to the velocities
            if( W > W_max ) {
              const CCTK_REAL W_scale = W_max/W;
              vxL = (vxL + betaxL)*W_scale - betaxL;
              vyL = (vyL + betayL)*W_scale - betayL;
              vzL = (vzL + betazL)*W_scale - betazL;
              W   = W_max;
            }

            // Step 3.f: Compute u^{mu} using:
            //  - W = alpha u^{0}     => u^{0} = W / alpha
            //  - v^{i} = u^{i}/u^{0} => u^{i} = v^{i}u^{0}
            const CCTK_REAL u0L = W / alpL;
            const CCTK_REAL uxL = vxL * u0L;
            const CCTK_REAL uyL = vyL * u0L;
            const CCTK_REAL uzL = vzL * u0L;

            // Step 3.g: Compute u_{mu} = g_{mu nu}u^{nu}
            // Step 3.g.i: Set g_{mu nu}
            // Step 3.g.i.1: Set gamma_{ij}
            CCTK_REAL gammaDD[3][3];
            gammaDD[0][0] = gxxL;
            gammaDD[0][1] = gammaDD[1][0] = gxyL;
            gammaDD[0][2] = gammaDD[2][0] = gxzL;
            gammaDD[1][1] = gyyL;
            gammaDD[1][2] = gammaDD[2][1] = gyzL;
            gammaDD[2][2] = gzzL;

            // Step 3.g.i.2: Compute beta_{i}
            CCTK_REAL betaU[3] = {betaxL,betayL,betazL};
            CCTK_REAL betaD[3] = {0.0,0.0,0.0};
            for(int ii=0;ii<3;ii++)
              for(int jj=0;jj<3;jj++)
                betaD[ii] += gammaDD[ii][jj] * betaU[jj];

            // Step 3.g.i.3: Compute beta^{2} = beta_{i}beta^{i}
            CCTK_REAL betasqr = 0.0;
            for(int ii=0;ii<3;ii++) betasqr += betaD[ii]*betaU[ii];

            // Step 3.g.i.4: Set g_{mu nu}
            CCTK_REAL g4DD[4][4];
            g4DD[0][0] = -alpL*alpL + betasqr;
            for(int ii=0;ii<3;ii++) {
              g4DD[0][ii+1] = g4DD[ii+1][0] = betaD[ii];
              for(int jj=ii;jj<3;jj++) {
                g4DD[ii+1][jj+1] = g4DD[jj+1][ii+1] = gammaDD[ii][jj];
              }
            }

            // Step 3.g.ii: Compute u_{mu} = g_{mu nu}u^{nu}
            CCTK_REAL u4U[4] = {u0L,uxL,uyL,uzL};
            CCTK_REAL u4D[4] = {0.0,0.0,0.0,0.0};
            for(int mu=0;mu<4;mu++)
              for(int nu=0;nu<4;nu++)
                u4D[mu] += g4DD[mu][nu] * u4U[nu];

            CCTK_REAL usqr = 0.0;
            for(int mu=0;mu<4;mu++) usqr += u4D[mu]*u4U[mu];
            if( fabs(usqr + 1) > 1e-13 ) {
              CCTK_VWARN(CCTK_WARN_ALERT,"Found u^2 != -1: usqr = %.15e",usqr);
            }

            // Step 3.h: Compute R, Q, and the neutrino opacities
            neutrino_opacities kappaL;
            CCTK_REAL R_sourceL, Q_sourceL;
            NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms(grhayl_eos,
                                                                          rhoL, Y_eL, temperatureL,
                                                                          &tauL, &kappaL, &R_sourceL, &Q_sourceL);

            if( robust_isnan(R_sourceL*Q_sourceL*kappaL.nue[0]*kappaL.nue[1]*kappaL.anue[0]*kappaL.anue[1]*kappaL.nux[0]*kappaL.nux[1]*
                             Ye_star_rhs[index]*tau_rhs[index]*st_x_rhs[index]*st_y_rhs[index]*st_z_rhs[index]) ) {
              CCTK_VINFO("****************************");
              CCTK_VINFO("NAN found:");
              CCTK_VINFO("rho Ye T: %e %e %e",rhoL,Y_eL,temperatureL);
              CCTK_VINFO("vx vy vz: %e %e %e",vxL,vyL,vzL);
              CCTK_VINFO("u^{mu}  : %e %e %e %e",u0L,uxL,uyL,uzL);
              CCTK_VINFO("alp beta: %e , %e %e %e",alpL,betaxL,betayL,betazL);
              CCTK_VINFO("R, Q | kappas: %e %e | %e %e , %e %e , %e %e",
                         R_sourceL,Q_sourceL,kappaL.nue[0],kappaL.nue[1],kappaL.anue[0],kappaL.anue[1],kappaL.nux[0],kappaL.nux[1]);
              CCTK_VINFO("rhss: %e %e %e %e %e",Ye_star_rhs[index],tau_rhs[index],st_x_rhs[index],st_y_rhs[index],st_z_rhs[index]);
              CCTK_VINFO("****************************");
              nan_found++;
            }

            // Step 3.i: Compute MHD right-hand sides
            const CCTK_REAL sqrtmgL      = alpL * psi6L;
            const CCTK_REAL sqrtmgR      = sqrtmgL * R_sourceL;
            const CCTK_REAL sqrtmgQ      = sqrtmgL * Q_sourceL;
            const CCTK_REAL Ye_star_rhsL = sqrtmgR;
            const CCTK_REAL tau_rhsL     = sqrtmgQ * u4U[0];
            const CCTK_REAL st_x_rhsL    = sqrtmgQ * u4D[1];
            const CCTK_REAL st_y_rhsL    = sqrtmgQ * u4D[2];
            const CCTK_REAL st_z_rhsL    = sqrtmgQ * u4D[3];

            // Step 3.j: Write to main memory
            kappa_0_nue [index]  = kappaL.nue [0];
            kappa_1_nue [index]  = kappaL.nue [1];
            kappa_0_anue[index]  = kappaL.anue[0];
            kappa_1_anue[index]  = kappaL.anue[1];
            kappa_0_nux [index]  = kappaL.nux [0];
            kappa_1_nux [index]  = kappaL.nux [1];

            // Step 3.k: Update right-hand sides only in the grid interior
            Ye_star_rhs [index] += Ye_star_rhsL;
            tau_rhs     [index] += tau_rhsL;
            st_x_rhs    [index] += st_x_rhsL;
            st_y_rhs    [index] += st_y_rhsL;
            st_z_rhs    [index] += st_z_rhsL;

            Ye_star_rhs_avg += Ye_star_rhsL;
            tau_rhs_avg     += tau_rhsL;
            st_x_rhs_avg    += st_x_rhsL;
            st_y_rhs_avg    += st_y_rhsL;
            st_z_rhs_avg    += st_z_rhsL;

            num_points++;
          }
        }
      }
    }
  }
  CCTK_REAL inv_numpts = num_points > 0 ? 1.0/((CCTK_REAL)num_points) : 1.0;
  if(verbosity_level>0) {
    CCTK_VINFO("***** Iter. # %d, Lev: %d, Averages -- Ye_rhs: %e | tau_rhs: %e | st_i_rhs: %e,%e,%e *****",cctk_iteration,GetRefinementLevel(cctkGH),
               Ye_star_rhs_avg*inv_numpts,
               tau_rhs_avg*inv_numpts,
               st_x_rhs_avg*inv_numpts,
               st_y_rhs_avg*inv_numpts,
               st_z_rhs_avg*inv_numpts);
    if(verbosity_level>1) CCTK_INFO("Finished NRPyLeakageET_compute_opacities_and_add_source_terms_to_MHD_rhss");
  }
  if( nan_found ) CCTK_ERROR("NAN Found. See error messages above. ABORTING!");
}
