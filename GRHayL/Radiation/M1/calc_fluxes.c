#include "ghl.h"
#include "ghl_m1.h"
#include "ghl_radiation.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

double minmod2(double rl, double rp, double th) {
    return fmin(1.0, fmin(th*rl, th*rp));
}

#define NGROUPS 1
#define NSPECIES 1

#define NSTENCIL_CONS 5
#define NSTENCIL_FLUX 3
#define NSTENCIL_CMAX 3 
#define NCONS_VAR_SIZE 5 // rN, rFx, rFy, rFz, rE

#define STENCIL_CONS_SIZE NGROUPS * NSPECIES * NSTENCIL_CONS
#define CONS_SIZE NGROUPS * NSPECIES * NSTENCIL_CONS * NCONS_VAR_SIZE
#define FLUX_SIZE NGROUPS * NSPECIES * NSTENCIL_FLUX * NCONS_VAR_SIZE
#define CMAX_SIZE NGROUPS * NSPECIES * NSTENCIL_FLUX

// Add IDX macros after testing
#define IDX_CONS(ig, k, iv) ((ig)*NSTENCIL_CONS*NCONS_VAR_SIZE + (k)*NCONS_VAR_SIZE + (iv))
#define IDX_FLUX(ig, k, iv) ((ig)*NSTENCIL_FLUX*NCONS_VAR_SIZE + (k)*NCONS_VAR_SIZE + (iv))
#define IDX_CMAX(ig, k)     ((ig)*NSTENCIL_CMAX + (k))

// #define GFINDEX1D(__k, ig, iv) \
//  ((iv) + (ig)*5 + (__k)*(5*NGROUPS*NSPECIES))


// This file translates from thc_M1_calc_fluxes.cc
// The goal of this function is to calculate the flux term d/dx^i F_i and store it in rhs
// Currently only implementing this for grid-based codes
int ghl_m1_calc_flux(
    double timestep,
    const ghl_metric_quantities stencil_metric[NSTENCIL_FLUX],
    const ghl_ADM_aux_quantities stencil_adm_aux[NSTENCIL_FLUX],
    const ghl_primitive_quantities stencil_prims[NSTENCIL_FLUX],
    const ghl_m1_thc_params *thc_params,
    const double stencil_rE[STENCIL_CONS_SIZE],
    const ghl_radiation_flux_vector stencil_rF4[STENCIL_CONS_SIZE],
    const ghl_radiation_pressure_tensor stencil_rP4[STENCIL_CONS_SIZE],
    const double stencil_rN[STENCIL_CONS_SIZE],
    const double delta[3], // [dx,dy,dz]
    const double abs_1[2],
    const double scat_1[2],
    const int dir,
    ghl_m1_conservative_quantities *rhs_flux // output
    ){

    double rhs[NGROUPS*NSPECIES*5] = {0};

    // When parallelized, these have to be on the heap
    double stencil_cons[CONS_SIZE]   = {0};
    double stencil_flux[FLUX_SIZE]   = {0};
    double stencil_cmax[CMAX_SIZE]   = {0};


    double flux_jm[NGROUPS*NSPECIES*5] = {0};
    double flux_jp[NGROUPS*NSPECIES*5] = {0};

    // Grid data
    double idelta[3] = {0};
    for (int i = 0; i < 3; i++) {
        idelta[i] = 1.0/delta[i];
    }

    // TODO: Check for inf or NaN 
    // ----------------------------------------------
    // 1st pass compute the fluxes
    // . this could be made into a separate Cactus scheduled
    //   function to avoid some duplicated calculations
    //  |  .  |  .  |  .  |  .  |  .  |
    //    k-2   k-1    k    k+1   k+2 
    //    umm   um     u    up    upp
    //          fm     f    fp
    //            fm1/2 fp1/2
    // ----------------------------------------------

    
    for (int ig = 0; ig < NGROUPS*NSPECIES; ++ig) {
        //      |  .  |  .  |  .  |  .  |  .  |
        // __k     0     1     2     3     4
        // i4D     0     1     2     3     4
        for (int __k = 0; __k < NSTENCIL_CONS; __k++){
            int i4D = ig*NSTENCIL_CONS + __k;
            ghl_radiation_flux_vector cur_rF4 = stencil_rF4[i4D];
            stencil_cons[IDX_CONS(ig, __k, 0)] =  stencil_rN[i4D];
            stencil_cons[IDX_CONS(ig, __k, 1)] = cur_rF4.D[1];
            stencil_cons[IDX_CONS(ig, __k, 2)] = cur_rF4.D[2];
            stencil_cons[IDX_CONS(ig, __k, 3)] = cur_rF4.D[3];
            stencil_cons[IDX_CONS(ig, __k, 4)] =  stencil_rE[i4D];
        } // k=4 stencil for cons

        //      |  .  |  .  |  .  |  .  |  .  |
        // __k           0     1     2      
        // i4D     0     1     2     3     4
        for (int __k = 0; __k < NSTENCIL_FLUX; __k++){
            int i4D = ig*NSTENCIL_CONS + __k + 1; // for cons total size is 4, there is a gap of 1
            double cur_rE = stencil_rE[i4D];
            ghl_radiation_flux_vector cur_rF4 = stencil_rF4[i4D];
            ghl_radiation_pressure_tensor cur_rP4 = stencil_rP4[i4D];

            double u4U[4] = { stencil_prims[__k].u0, \
                stencil_prims[__k].u0 * stencil_prims[__k].vU[0], \
                stencil_prims[__k].u0 * stencil_prims[__k].vU[1], \
                stencil_prims[__k].u0 * stencil_prims[__k].vU[2] };
            double v4U[4] = { 0, stencil_prims[__k].vU[0], stencil_prims[__k].vU[1], stencil_prims[__k].vU[2] };

            double u4D[4] = { 0, 0, 0, 0 };
            double v4D[4] = { 0, 0, 0, 0 };
            for(int a = 0; a < 4; a++) {
                for(int b = 0; b < 4; b++) {
                    u4D[a] += stencil_adm_aux[__k].g4DD[a][b] * u4U[b];
                    v4D[a] += stencil_adm_aux[__k].g4DD[a][b] * v4U[b];
                }
            }
            double n4D[4] = { -stencil_metric[__k].lapse, 0, 0, 0 };
            double n4U[4]
            = { stencil_metric[__k].lapseinv, -stencil_metric[__k].betaU[0] * stencil_metric[__k].lapseinv,
                                            -stencil_metric[__k].betaU[1] * stencil_metric[__k].lapseinv,
                                            -stencil_metric[__k].betaU[2] * stencil_metric[__k].lapseinv };
            double W = 0;
            for(int a = 0; a < 4; a++) {
                W += -u4D[a] * n4U[a];
            }
            ghl_radiation_metric_tensor proj4;
            for(int a = 0; a < 4; a++) {
                for(int b = 0; b < 4; b++) {
                    proj4.UD[a][b] = (int)(a == b) + u4U[a] * u4D[b];
                }
            }
            ghl_stress_energy rT4DD = {0};
            assemble_rT_lab_frame(n4D, cur_rE, &cur_rF4, &cur_rP4, &rT4DD);
            double J = calc_J_from_rT(u4U, &proj4, &rT4DD);
            ghl_radiation_flux_vector H4;
            calc_H4D_from_rT(u4U, &proj4, &rT4DD, &H4);


            double Gamma = compute_Gamma(W, &v4D, J, cur_rE,
                thc_params->rad_E_floor, thc_params->rad_eps,
                &cur_rF4);

            ghl_radiation_con_flux_vector fnu4 = {0};
            assemble_fnu(&stencil_adm_aux[__k], &u4U, J, &H4, &fnu4);

            // Note that nnu is densitized here
            double nnu = stencil_rN[__k]/Gamma;
            // calc_F_flux uses F4.U
            get_FU(&stencil_adm_aux[__k], &cur_rF4);
            stencil_flux[IDX_FLUX(ig, __k, 0)] = stencil_metric[__k].lapse * nnu * fnu4.U[dir+1];
            stencil_flux[IDX_FLUX(ig, __k, 1)] = calc_F_flux(&stencil_metric[__k], &cur_rF4, &cur_rP4, dir+1, 1);
            stencil_flux[IDX_FLUX(ig, __k, 2)] = calc_F_flux(&stencil_metric[__k], &cur_rF4, &cur_rP4, dir+1, 2);
            stencil_flux[IDX_FLUX(ig, __k, 3)] = calc_F_flux(&stencil_metric[__k], &cur_rF4, &cur_rP4, dir+1, 3);
            stencil_flux[IDX_FLUX(ig, __k, 4)] = calc_E_flux(&stencil_metric[__k], cur_rE, &cur_rF4, dir+1);
            

            // Speed of light -- note that gammaUU has NDIM not 3 from THC but 4 here
            double clam_1 =  stencil_metric[__k].lapse*sqrt(stencil_metric[__k].gammaUU[dir+1][dir+1]) - stencil_metric[__k].betaU[dir+1];
            double clam_2 = -stencil_metric[__k].lapse*sqrt(stencil_metric[__k].gammaUU[dir+1][dir+1]) - stencil_metric[__k].betaU[dir+1];
            // In some cases the eigenvalues in the thin limit
            // overestimate the actual characteristic speed and can
            // become larger than c
            double clight = fmax(fabs(clam_1), fabs(clam_2));
            stencil_cmax[IDX_CMAX(ig, __k)] = clight;
                                                 // = std::min(clight, cM1);
        } // k=2 stencil for flux and cmax

        // ----------------------------------------------
        // 2nd pass update the RHS
        for (int iv = 0; iv < NSTENCIL_CONS; ++iv) {
            double const ujmm =  stencil_cons[IDX_CONS(ig, 0, iv)];
            double const ujm  =  stencil_cons[IDX_CONS(ig, 1, iv)];
            double const uj   =  stencil_cons[IDX_CONS(ig, 2, iv)];
            double const ujp  =  stencil_cons[IDX_CONS(ig, 3, iv)];
            double const ujpp =  stencil_cons[IDX_CONS(ig, 4, iv)];

            double const fjm  =  stencil_flux[IDX_FLUX(ig, 0, iv)];
            double const fj   =  stencil_flux[IDX_FLUX(ig, 1, iv)];
            double const fjp  =  stencil_flux[IDX_FLUX(ig, 2, iv)];

            double const ccm  =  stencil_cmax[IDX_CMAX(ig, 0)];
            double const cc   =  stencil_cmax[IDX_CMAX(ig, 1)];
            double const ccp  =  stencil_cmax[IDX_CMAX(ig, 2)];
            double const cmx =  fmax(cc, ccp);

            double const dupp = ujpp - ujp;
            double const dupc = ujp - uj;
            double const ducm = uj - ujm;
            double const dumm = uj - ujm;

            // Remove dissipation at high Peclet numbers
            double kapa = 0.5*(abs_1[1] + abs_1[0] +
                    scat_1[1] + scat_1[0]);
            double A = 1.0;
            if (kapa*delta[dir] > 1.0) {
                A = fmin(1.0, 1.0/(delta[dir]*kapa));
                A = fmax(A, thc_params->mindiss);
            }

            // Instead of rotating pointers like THC, we just calculate both flux_jm and flux_jp here
            // flux_jp
            bool sawtooth_p = false;
            double phi_p = 0;
            if (dupc*dupp > 0 && dupc*ducm > 0) {
                phi_p = minmod2(dumm/ducm, dupc/dupp, thc_params->minmod_theta);
            }
            else if (dupc*dupp < 0 && dupc*ducm < 0) {
                sawtooth_p = true;
            }
            //assert(isfinite(phi));
            double const flux_p_low =
            0.5*(fj + fjp - cmx*(ujp - uj));
            double const flux_p_high = 0.5*(fj + fjp);

            double flux_p_num = flux_p_high -
                (sawtooth_p ? 1 : A)*(1 - phi_p)*(flux_p_high - flux_p_low);

            // flux_jm
            bool sawtooth_m = false;
            double phi_m = 0;
            if (ducm*dupc > 0 && ducm*dumm > 0) {
                phi_m = minmod2(dumm/ducm, dupc/dupp, thc_params->minmod_theta);
            }
            else if (ducm*dupc < 0 && ducm*dumm < 0) {
                sawtooth_m = true;
            }
            //assert(isfinite(phi));
            double const flux_m_low =
            0.5*(fjm + fj - cmx*(uj - ujm));
            double const flux_m_high = 0.5*(fjm + fj);

            double flux_m_num = flux_m_high -
                (sawtooth_m ? 1 : A)*(1 - phi_m)*(flux_m_high - flux_m_low);

            int i4D_rhs = ig*5 + iv;
            flux_jp[i4D_rhs] = flux_p_num;
            flux_jm[i4D_rhs] = flux_m_num;

            rhs[i4D_rhs] += idelta[dir]*(
                flux_jm[i4D_rhs] -
                flux_jp[i4D_rhs]);
        } // for iv

        rhs_flux->rN[ig]  = rhs[ig*5 + 0];
        rhs_flux->rF4[ig].D[0] = 0.0; // No source for rF4[ig].D[0] in fluxes
        rhs_flux->rF4[ig].D[1] = rhs[ig*5 + 1];
        rhs_flux->rF4[ig].D[2] = rhs[ig*5 + 2];
        rhs_flux->rF4[ig].D[3] = rhs[ig*5 + 3];
        rhs_flux->rE[ig]  = rhs[ig*5 + 4];
    } // for ig
    return 1;
} // ghl_m1_calc_flux