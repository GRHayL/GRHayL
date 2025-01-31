#include "ghl.h"
#include "ghl_radiation.h"

// This file translates from thc_M1_calc_update.cc

void ghl_M1_update(
    double cdt,
    ghl_neutrino_optical_depths *restrict tau,
    ghl_neutrino_opacities *restrict kappa,
    ghl_metric_quantities *metric,
    ghl_ADM_aux_quantities *adm_aux,
    ghl_radiation_metric_tensor *proj4,
    ghl_primitive_quantities *prims
    ){
    // Disable GSL error handler
    gsl_error_handler_t * gsl_err = gsl_set_error_handler_off();
    // Set the closure type
    ghl_m1_closure_t ghl_m1_closure_type = Minerbo; //TODO: is this right?

    /////////// THC params: ////////////
    // source_limiter: Limiter value (0: sources are disabled, 1: sources are limited 
    double source_limiter = -1;
    // source_Ye_min/max: Minimum/Maximum allowed Ye for the matter"
    double source_Ye_min = 0.0;
    double source_Ye_max = 0.6;
    double rad_E_floor = 1.0e-15;
    double rad_eps =   1.0e-15;
    double source_therm_limit = -1.0;
    int ngroups = 1;
    int nspecies = 3;
    ///////////////////////////////////


    // Steps
    // 1. F^m   = F^k + dt/2 [ A[F^k] + S[F^m]   ]
    // 2. F^k+1 = F^k + dt   [ A[F^m] + S[F^k+1] ]
    // At each step we solve an implicit problem in the form
    //    F = F^* + cdt S[F]
    // Where F^* = F^k + cdt A
    
    //int const siz = UTILS_GFSIZE(cctkGH);
    double * sconx = &scon[0*siz];
    double * scony = &scon[1*siz];
    double * sconz = &scon[2*siz];

    double const mb = AverageBaryonMass();

    // begin parallel
    // do these come from prims?
    double dens[ngroups*nspecies];
    double Y_e[ngroups*nspecies];
    double nueave[ngroups*nspecies];
    double chi[ngroups*nspecies];

    double fidu_w_lorentz = 1.0; //TODO
    double eta_1 = 1.0; //TODO
    double abs_1[ngroups*nspecies];
    double volform = 1.0; //TODO
    double rad_N_floor = 1.0; //TODO
    double scat_1[ngroups*nspecies];
    double rE_p[ngroups*nspecies];
    double rFx_p[ngroups*nspecies];
    double rFy_p[ngroups*nspecies];
    double rFz_p[ngroups*nspecies];
    double rN_p[ngroups*nspecies];
    double rE_rhs[ngroups*nspecies];
    double rFx_rhs[ngroups*nspecies];
    double rFy_rhs[ngroups*nspecies];
    double rFz_rhs[ngroups*nspecies];
    double rN_rhs[ngroups*nspecies];
    double E[ngroups*nspecies];
    double F4[ngroups*nspecies][4];
    double H4_star[ngroups*nspecies][4];
    double g_uu[4][4];
    double u4U[4];
    
    double rT4DD[4][4];
    double n4D[4];
    double n4U[4];

    //
    // Source RHS are stored here
    double DrE[ngroups*nspecies];
    double DrFx[ngroups*nspecies];
    double DrFy[ngroups*nspecies];
    double DrFz[ngroups*nspecies];
    double DrN[ngroups*nspecies];
    double DDxp[ngroups*nspecies];

    // Step 1 -- compute the sources
    for (int ig = 0; ig < ngroups*nspecies; ++ig) {
        //
        // Advect radiation
        double E_star = rE_p[ig] + cdt*rE_rhs[ig];
        double F4_star[4];
        F4_star[1] = rFx_p[ig] + cdt*rFx_rhs[ig];
        F4_star[2] = rFy_p[ig] + cdt*rFy_rhs[ig];
        F4_star[3] = rFz_p[ig] + cdt*rFz_rhs[ig];
        apply_floor(adm_aux, E_star, &F4_star, rad_E_floor, rad_eps);
        double N_star = max(rN_p[ig] + cdt*rN_rhs[ig], rad_N_floor);
        double E_new;

        //
        // Compute quantities in the fluid frame
        ghl_radiation_pressure_tensor P4;
        m1_root_params closure_params;
        closure_params.metric   = &metric;
        closure_params.adm_aux  = &adm_aux;
        closure_params.prims    = &prims;
        closure_params.E        = E;
        closure_params.F4       = &F4;
        closure_params.P4       = &P4;
        ghl_radiation_rootSolve_closure(closure_params);
        // Need to calculate if P4 is inf

        ghl_stress_energy *rT4DD;
        assemble_rT_lab_frame(n4D, E_star, F4_star, P4, &rT4DD);

        double const Jstar = calc_J_from_rT(u4U, proj4, rT4DD);

        ghl_radiation_flux_vector H4_star;
        calc_H_from_rT(u4U, proj4, rT4DD, &H4_star);

        //
        // Estimate interaction with matter
        double const dtau = cdt/fidu_w_lorentz;
        double J_new = (Jstar + dtau*eta_1*volform)/(1 + dtau*abs_1[ig]);

        // Only three components of H^a are independent H^0 is found by
        // requiring H^a u_a = 0
        double const khat = (abs_1[ig] + scat_1[ig]);
        ghl_radiation_flux_vector * Hnew4;
        for (int a = 1; a < 4; ++a) {
            Hnew4->D[a] = H4_star->D[a]/(1 + dtau*khat);
        }
        Hnew4->D[0] = 0.0;
        for (int a = 1; a < 4; ++a) {
            Hnew4->D[0] -= Hnew4->D[a]*(u4U[a]/u4U[0]);
        }

        //
        // Update Tmunu
        double const H2 = ghl_compute_vec2_from_vec4D(adm_aux->g4UU, Hnew_d, Hnew_d);
        double const dthick = 3.*(1. - chi[ig])/2.;
        double const dthin = 1. - dthick;

        for(int a = 0; a < 4; ++a) {
            for(int b = 0; b < 4; ++b) {
                rT4DD[a][b] = J_new*u4U[a]*u4U[b] + Hnew4->D[a]*u4U[b] + Hnew4->D[b]*u4U[a] +
                    dthin*J_new*(Hnew4->D[a]*Hnew4->D[b]*(H2 > 0 ? 1/H2 : 0)) +
                    dthick*J_new*(adm_aux->g4DD[a][b] + u4U[a]*u4U[b])/3;
            }
        }

        //
        // Boost back to the lab frame
        E_new = calc_J_from_rT(rT4DD, n4U);
        calc_H4D_from_rT(u4U, &proj4, &rT4DD, &H4_star);
        apply_floor(adm_aux, E_new, &F4_new, rad_E_floor, rad_eps);
        //
        // Compute interaction with matter
        (void)source_update(); //TODO!
        apply_floor(adm_aux, E_new, &F4_new, rad_E_floor, rad_eps);

        //
        // Update closure
        ghl_radiation_apply_closure(metric, adm_aux, prims, E_new, F4_new, chi[ig], &P4);

        //
        // Compute new radiation energy density in the fluid frame
        ghl_stress_energy T4;
        assemble_rT(n4D, E_new, &F4_new, &P4, &T4);
        J_new = calc_J_from_rT(T4, u4U);

        //
        // Compute changes in radiation energy and momentum
        DrE[ig]  = E_new - E_star;
        DrFx[ig] = F4_new->D[1] - F4_star->D[1];
        DrFy[ig] = F4_new->D[2] - F4_star->D[2];
        DrFz[ig] = F4_new->D[3] - F4_star->D[3];

        //
        // Compute updated Gamma
        double const Gamma = compute_Gamma(fidu_w_lorentz, &v4U, J_new, E_new, rad_E_floor, rad_eps, F4_new);

        //
        // N^k+1 = N^* + dt ( eta - abs N^k+1 )
        if (source_therm_limit < 0 || dt*abs_0[ig] < source_therm_limit) {
            DrN[ig] = (N_star + dt*metric->lapse*volform*eta_0[ig])/
                        (1 + dt*metric->lapse*abs_0[ig]/Gamma) - N_star;
        }
        //
        // The neutrino number density is updated assuming the neutrino
        // average energies are those of the equilibrium
        else {
            DrN[ig] = (nueave[ig] > 0 ? Gamma*J_new/nueave[ig] - N_star : 0.0);
        }

        //
        // Fluid lepton sources
        DDxp[ig] = -mb*(DrN[ig]*(ig == 0) - DrN[ig]*(ig == 1));
    } // for (int ig ...)
    //
    // Step 2 -- limit the sources
    double theta = 1.0;
    

    if (source_limiter >= 0) {
        theta = 1.0;
        double DTau_sum = 0.0;
        double DDxp_sum = 0.0;
        for (int ig = 0; ig < ngroups*nspecies; ++ig) {
            double E_star = rE_p[ig] + cdt*rE_rhs[ig];
            if (DrE[ig] < 0) {
                theta = min(-source_limiter*max(E_star, 0.0)/DrE[ig], theta);
            }
            DTau_sum -= DrE[ig];

            double N_star = rN_p[ig] + cdt*rN_rhs[ig];
            if (DrN[ig] < 0) {
                theta = min(-source_limiter*max(N_star, 0.0)/DrN[ig], theta);
            }
            DDxp_sum += DDxp[ig];
        }
        double const DYe = DDxp_sum/dens;
        if (DTau_sum < 0) {
            theta = min(-source_limiter*max(tau, 0.0)/DTau_sum, theta);
        }
        if (DYe > 0) {
            theta = min(source_limiter*max(source_Ye_max - Y_e, 0.0)/DYe, theta);
        }
        else if (DYe < 0) {
            theta = min(source_limiter*min(source_Ye_min - Y_e, 0.0)/DYe, theta);
        }
        theta = max(0.0, theta);
    }

    //
    // Step 3 -- update fields
    for (int ig = 0; ig < ngroups*nspecies; ++ig) {

        //
        // Update radiation quantities
        double E =  rE_p[ig] + cdt*rE_rhs[ig]  + theta*DrE[ig];
        F4->D[1]      = rFx_p[ig] + cdt*rFx_rhs[ig] + theta*DrFx[ig];
        F4->D[2]      = rFy_p[ig] + cdt*rFy_rhs[ig] + theta*DrFy[ig];
        F4->D[3]      = rFz_p[ig] + cdt*rFz_rhs[ig] + theta*DrFz[ig];
        apply_floor(adm_aux, &E, &F4, rad_E_floor, rad_eps);

        double N =  rN_p[ig] + cdt*rN_rhs[ig]  + theta*DrN[ig];
        N = max(N, rad_N_floor);

        //
        // Compute back reaction on the fluid
        // NOTE: fluid backreaction is only needed at the last substep
        if (backreact && 0 == *TimeIntegratorStage) {
            assert (ngroups == 1);
            assert (nspecies == 3);

            sconx  -= theta*DrFx[ig];
            scony  -= theta*DrFy[ig];
            sconz  -= theta*DrFz[ig];
            tau    -= theta*DrE[ig];
            densxp += theta*DDxp[ig];
            densxn -= theta*DDxp[ig];

            netabs  += theta*DDxp[ig];
            netheat -= theta*DrE[ig];
        }

        //
        // Save updated results into grid functions
        rE[ig]  = E;
        unpack_F_d(F_d, &rFx[ig], &rFy[ig], &rFz[ig]);
        rN[ig] = N;
    }
    gsl_root_fsolver_free(gsl_solver_1d);
    gsl_multiroot_fdfsolver_free(gsl_solver_nd);
    // End parallel 

    // Restore GSL error handler
    gsl_set_error_handler(gsl_err);
} // ghl_M1_update