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

  const double poison = 0.0/0.0;
  double dummy1, dummy2, dummy3;

  const int imax = cctk_lsh[0];
  const int jmax = cctk_lsh[1];
  const int kmax = cctk_lsh[2];

// We use rho and press from HydroBase directly with no need to convert
#pragma omp parallel for
  for(int k=0; k<kmax; k++) {
    for(int j=0; j<jmax; j++) {
      for(int i=0; i<imax; i++) {
        const int index=CCTK_GFINDEX3D(cctkGH,i,j,k);
        const int ind0=CCTK_GFINDEX4D(cctkGH,i,j,k,0);
        const int ind1=CCTK_GFINDEX4D(cctkGH,i,j,k,1);
        const int ind2=CCTK_GFINDEX4D(cctkGH,i,j,k,2);

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

  // Next compute B & B_stagger from A_i. Note that this routine also depends on
  //   the psi_bssn[] gridfunction being set to exp(phi).

  const double dxi = 1.0/CCTK_DELTA_SPACE(0);
  const double dyi = 1.0/CCTK_DELTA_SPACE(1);
  const double dzi = 1.0/CCTK_DELTA_SPACE(2);

#pragma omp parallel for
  for(int k=0; k<kmax; k++) {
    for(int j=0; j<jmax; j++) {
      for(int i=0; i<imax; i++) {
        // Look Mom, no if() statements!
        const int shiftedim1 = (i-1)*(i!=0); // This way, i=0 yields shiftedim1=0 and shiftedi=1, used below for our COPY boundary condition.
        const int shiftedi   = shiftedim1+1;

        const int shiftedjm1 = (j-1)*(j!=0);
        const int shiftedj   = shiftedjm1+1;

        const int shiftedkm1 = (k-1)*(k!=0);
        const int shiftedk   = shiftedkm1+1;

        int index,indexim1,indexjm1,indexkm1;

        const int actual_index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        const double Psi = psi_bssn[actual_index];
        const double Psim3 = 1.0/(Psi*Psi*Psi);

        // For the lower boundaries, the following applies a "copy"
        //    boundary condition on Bi_stagger where needed.
        //    E.g., Bx_stagger(i,jmin,k) = Bx_stagger(i,jmin+1,k)
        //    We find the copy BC works better than extrapolation.
        // For the upper boundaries, we do the following copy:
        //    E.g., Psi(imax+1,j,k)=Psi(imax,j,k)
        /**************/
        /* Bx_stagger */
        /**************/

        index    = CCTK_GFINDEX3D(cctkGH,i,shiftedj,shiftedk);
        indexjm1 = CCTK_GFINDEX3D(cctkGH,i,shiftedjm1,shiftedk);
        indexkm1 = CCTK_GFINDEX3D(cctkGH,i,shiftedj,shiftedkm1);
        // Set Bx_stagger = \partial_y A_z - partial_z A_y
        // "Grid" Ax(i,j,k) is actually Ax(i,j+1/2,k+1/2)
        // "Grid" Ay(i,j,k) is actually Ay(i+1/2,j,k+1/2)
        // "Grid" Az(i,j,k) is actually Ay(i+1/2,j+1/2,k)
        // Therefore, the 2nd order derivative \partial_z A_y at (i+1/2,j,k) is:
        //          ["Grid" Ay(i,j,k) - "Grid" Ay(i,j,k-1)]/dZ
        Bx_stagger[actual_index] = (Az[index]-Az[indexjm1])*dyi - (Ay[index]-Ay[indexkm1])*dzi;

        // Now multiply Bx and Bx_stagger by 1/sqrt(gamma(i+1/2,j,k)]) = 1/sqrt(1/2 [gamma + gamma_ip1]) = exp(-6 x 1/2 [phi + phi_ip1] )
        const int imax_minus_i = (cctk_lsh[0]-1)-i;
        const int indexip1jk = CCTK_GFINDEX3D(cctkGH,i + ( (imax_minus_i > 0) - (0 > imax_minus_i) ),j,k);
        const double Psi_ip1 = psi_bssn[indexip1jk];
        Bx_stagger[actual_index] *= Psim3/(Psi_ip1*Psi_ip1*Psi_ip1);

        /**************/
        /* By_stagger */
        /**************/

        index    = CCTK_GFINDEX3D(cctkGH,shiftedi,j,shiftedk);
        indexim1 = CCTK_GFINDEX3D(cctkGH,shiftedim1,j,shiftedk);
        indexkm1 = CCTK_GFINDEX3D(cctkGH,shiftedi,j,shiftedkm1);
        // Set By_stagger = \partial_z A_x - \partial_x A_z
        By_stagger[actual_index] = (Ax[index]-Ax[indexkm1])*dzi - (Az[index]-Az[indexim1])*dxi;

        // Now multiply By and By_stagger by 1/sqrt(gamma(i,j+1/2,k)]) = 1/sqrt(1/2 [gamma + gamma_jp1]) = exp(-6 x 1/2 [phi + phi_jp1] )
        const int jmax_minus_j = (cctk_lsh[1]-1)-j;
        const int indexijp1k = CCTK_GFINDEX3D(cctkGH,i,j + ( (jmax_minus_j > 0) - (0 > jmax_minus_j) ),k);
        const double Psi_jp1 = psi_bssn[indexijp1k];
        By_stagger[actual_index] *= Psim3/(Psi_jp1*Psi_jp1*Psi_jp1);


        /**************/
        /* Bz_stagger */
        /**************/

        index    = CCTK_GFINDEX3D(cctkGH,shiftedi,shiftedj,k);
        indexim1 = CCTK_GFINDEX3D(cctkGH,shiftedim1,shiftedj,k);
        indexjm1 = CCTK_GFINDEX3D(cctkGH,shiftedi,shiftedjm1,k);
        // Set Bz_stagger = \partial_x A_y - \partial_y A_x
        Bz_stagger[actual_index] = (Ay[index]-Ay[indexim1])*dxi - (Ax[index]-Ax[indexjm1])*dyi;

        // Now multiply Bz_stagger by 1/sqrt(gamma(i,j,k+1/2)]) = 1/sqrt(1/2 [gamma + gamma_kp1]) = exp(-6 x 1/2 [phi + phi_kp1] )
        const int kmax_minus_k = (cctk_lsh[2]-1)-k;
        const int indexijkp1 = CCTK_GFINDEX3D(cctkGH,i,j,k + ( (kmax_minus_k > 0) - (0 > kmax_minus_k) ));
        const double Psi_kp1 = psi_bssn[indexijkp1];
        Bz_stagger[actual_index] *= Psim3/(Psi_kp1*Psi_kp1*Psi_kp1);
      }
    }
  }

#pragma omp parallel for
  for(int k=0; k<kmax; k++) {
    for(int j=0; j<jmax; j++) {
      for(int i=0; i<imax; i++) {
        // Look Mom, no if() statements!
        const int shiftedim1 = (i-1)*(i!=0); // This way, i=0 yields shiftedim1=0 and shiftedi=1, used below for our COPY boundary condition.
        const int shiftedi   = shiftedim1+1;

        const int shiftedjm1 = (j-1)*(j!=0);
        const int shiftedj   = shiftedjm1+1;

        const int shiftedkm1 = (k-1)*(k!=0);
        const int shiftedk   = shiftedkm1+1;

        int index,indexim1,indexjm1,indexkm1;

        const int actual_index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        // For the lower boundaries, the following applies a "copy"
        //    boundary condition on Bi and Bi_stagger where needed.
        //    E.g., Bx(imin,j,k) = Bx(imin+1,j,k)
        //    We find the copy BC works better than extrapolation.
        /******/
        /* Bx */
        /******/
        index = CCTK_GFINDEX3D(cctkGH,shiftedi,j,k);
        indexim1 = CCTK_GFINDEX3D(cctkGH,shiftedim1,j,k);
        // Set Bx = 0.5 ( Bx_stagger + Bx_stagger_im1 )
        // "Grid" Bx_stagger(i,j,k) is actually Bx_stagger(i+1/2,j,k)
        Bx_center[actual_index] = 0.5 * ( Bx_stagger[index] + Bx_stagger[indexim1] );

        /******/
        /* By */
        /******/
        index = CCTK_GFINDEX3D(cctkGH,i,shiftedj,k);
        indexjm1 = CCTK_GFINDEX3D(cctkGH,i,shiftedjm1,k);
        // Set By = 0.5 ( By_stagger + By_stagger_im1 )
        // "Grid" By_stagger(i,j,k) is actually By_stagger(i,j+1/2,k)
        By_center[actual_index] = 0.5 * ( By_stagger[index] + By_stagger[indexjm1] );

        /******/
        /* Bz */
        /******/
        index = CCTK_GFINDEX3D(cctkGH,i,j,shiftedk);
        indexkm1 = CCTK_GFINDEX3D(cctkGH,i,j,shiftedkm1);
        // Set Bz = 0.5 ( Bz_stagger + Bz_stagger_im1 )
        // "Grid" Bz_stagger(i,j,k) is actually Bz_stagger(i,j+1/2,k)
        Bz_center[actual_index] = 0.5 * ( Bz_stagger[index] + Bz_stagger[indexkm1] );
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
        grhayl_initialize_metric(
              alp[index],
              betax[index], betay[index], betaz[index],
              gxx[index], gxy[index], gxz[index],
              gyy[index], gyz[index], gzz[index],
              &ADM_metric);

        ADM_aux_quantities metric_aux;
        grhayl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

        primitive_quantities prims;
        grhayl_initialize_primitives(
              rho_b[index], pressure[index], eps[index],
              vx[index], vy[index], vz[index],
              Bx_center[index], By_center[index], Bz_center[index],
              poison, poison, poison,
              &prims);

        conservative_quantities cons;
        stress_energy Tmunu;
        int speed_limited = 0;
        //This applies inequality fixes on the conservatives
        grhayl_enforce_primitive_limits_and_compute_u0(
              grhayl_params, grhayl_eos, &ADM_metric,
              &metric_aux, &prims, &speed_limited);
        //This computes the conservatives and stress-energy tensor from the new primitives
        grhayl_compute_conservs_and_Tmunu(
              &ADM_metric, &metric_aux, &prims, &cons, &Tmunu);

        grhayl_return_primitives(
              &prims,
              &rho_b[index], &pressure[index], &eps[index],
              &vx[index], &vy[index], &vz[index],
              &Bx_center[index], &By_center[index], &Bz_center[index],
              &dummy1, &dummy2, &dummy3);

        grhayl_return_conservatives(
              &cons,
              &rho_star[index], &tau[index],
              &Stildex[index], &Stildey[index], &Stildez[index],
              &dummy1, &dummy2);

        if(update_Tmunu) {
          grhayl_return_stress_energy(
                &Tmunu, &eTtt[index], &eTtx[index],
                &eTty[index], &eTtz[index], &eTxx[index],
                &eTxy[index], &eTxz[index], &eTyy[index],
                &eTyz[index], &eTzz[index]);
        }
      }
    }
  }
}
