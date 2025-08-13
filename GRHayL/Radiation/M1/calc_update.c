#include "ghl.h"
#include "ghl_radiation.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

// This file translates from thc_M1_calc_update.cc

void ghl_M1_update(
    double timestep,
    const ghl_metric_quantities *metric,
    const ghl_metric_quantities *metric_derivs_x,
    const ghl_metric_quantities *metric_derivs_y,
    const ghl_metric_quantities *metric_derivs_z,
    const ghl_extrinsic_curvature *curv,
    const ghl_ADM_aux_quantities *adm_aux,
    const ghl_primitive_quantities *prims,
    const ghl_m1_thc_params *thc_params
    ){

    /////////// THC params: ////////////
    // source_limiter: Limiter value (0: sources are disabled, 1: sources are limited 
    double source_limiter = -1;
    // source_Ye_min/max: Minimum/Maximum allowed Ye for the matter"
    double source_Ye_min = 0.0;
    double source_Ye_max = 0.6;
    double rad_E_floor = 1.0e-15;
    double rad_eps =   1.0e-15;
    double rad_N_floor = 1.0e-10;
    double source_therm_limit = -1.0;
    #define NGROUPS 1
    #define NSPECIES 1
    ///////////////////////////////////


    // Public grid variables (in interface.ccl)
    // Evolved quantities in current time step
    double rE[NGROUPS*NSPECIES]    = {0};
    double rFx[NGROUPS*NSPECIES]   = {0};
    double rFy[NGROUPS*NSPECIES]   = {0};
    double rFz[NGROUPS*NSPECIES]   = {0};
    double rN[NGROUPS*NSPECIES]    = {0};

    double rPxx[NGROUPS*NSPECIES]  = {0};
    double rPxy[NGROUPS*NSPECIES]  = {0};
    double rPxz[NGROUPS*NSPECIES]  = {0};
    double rPyy[NGROUPS*NSPECIES]  = {0};
    double rPyz[NGROUPS*NSPECIES]  = {0};
    double rPzz[NGROUPS*NSPECIES]  = {0};

    // Evolved quantities in previous time step
    double rE_p[NGROUPS*NSPECIES]  = {0};
    double rFx_p[NGROUPS*NSPECIES] = {0};
    double rFy_p[NGROUPS*NSPECIES] = {0};
    double rFz_p[NGROUPS*NSPECIES] = {0};
    double rN_p[NGROUPS*NSPECIES]  = {0};

    // fluid frame
    double rJ[NGROUPS*NSPECIES]    = {0};
    double rHt[NGROUPS*NSPECIES]   = {0};
    double rHx[NGROUPS*NSPECIES]   = {0};
    double rHy[NGROUPS*NSPECIES]   = {0};
    double rHz[NGROUPS*NSPECIES]   = {0};

    double rE_rhs[NGROUPS*NSPECIES]   = {0};
    double rFx_rhs[NGROUPS*NSPECIES]  = {0};
    double rFy_rhs[NGROUPS*NSPECIES]  = {0};
    double rFz_rhs[NGROUPS*NSPECIES]  = {0};
    double rN_rhs[NGROUPS*NSPECIES]   = {0};

    double nueave[NSPECIES*NGROUPS]   = {0};
    double chi[NGROUPS*NSPECIES]      = {0};
        // opacities
    double abs_0[NGROUPS*NSPECIES]    = {0};
    double abs_1[NGROUPS*NSPECIES]    = {0};
    double eta_0[NGROUPS*NSPECIES]    = {0};
    double eta_1[NGROUPS*NSPECIES]    = {0};
    double scat_1[NGROUPS*NSPECIES]   = {0};

    // THC flag
    int M1_source_method = 0; // 0: explicit, 1: implicit, 2: boost
    int M1_source_method_explicit = 0;
    int M1_source_method_implicit = 1;
    int M1_source_method_boost    = 2;


    // THC_Core variables (in schedule.ccl)
    double volform = 1.0;
    double dens = 0.0;
    double tau = 0.0;
    double Y_e = 0.0;
    // Johnny: these are probably for back reaction
    // double* scon[1];
    // int const siz = UTILS_GFSIZE(cctkGH);
    // double * sconx = &scon[0*siz];
    // double * scony = &scon[1*siz];
    // double * sconz = &scon[2*siz];

    // Disable GSL error handler
    gsl_error_handler_t * gsl_err = gsl_set_error_handler_off();

    // Set the closure type
    ghl_m1_closure_t ghl_m1_closure_type = Minerbo;


    // Steps
    // 1. F^m   = F^k + dt/2 [ A[F^k] + S[F^m]   ]
    // 2. F^k+1 = F^k + dt   [ A[F^m] + S[F^k+1] ]
    // At each step we solve an implicit problem in the form
    //    F = F^* + cdt S[F]
    // Where F^* = F^k + cdt A

    double const mb = 0.0; //AverageBaryonMass(); //TODO: what exactly is this function?
    
    for(int timeIntegrator = 0; timeIntegrator < 2; ++timeIntegrator) {
        double cdt = (timeIntegrator == 0) ? 0.5*timestep : timestep;
        // begin parallel
        ghl_m1_closure_t closure = Minerbo;
        gsl_root_fsolver * gsl_solver_1d =
                gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
        gsl_multiroot_fdfsolver * gsl_solver_nd =
            gsl_multiroot_fdfsolver_alloc(gsl_multiroot_fdfsolver_hybridsj, 4);
        ghl_m1_powell_params powell_params = {0};
         
        ghl_stress_energy rT4DD;
        

        double E = 0.0;
        double E_star = 0.0;
        double E_new = 0.0;

        ghl_radiation_flux_vector F4 = {0};
        ghl_radiation_flux_vector F4_star = {0};
        ghl_radiation_flux_vector F4_new = {0};

        // Initialize P4
        ghl_radiation_pressure_tensor P4;
        for(int a = 0; a < 4; ++a) {
            for(int b = 0; b < 4; ++b) {
                P4.DD[a][b] = 0.0;
            }
        }

        // Source RHS are stored here
        double DrE[NGROUPS*NSPECIES];
        double DrFx[NGROUPS*NSPECIES];
        double DrFy[NGROUPS*NSPECIES];
        double DrFz[NGROUPS*NSPECIES];
        double DrN[NGROUPS*NSPECIES];
        double DDxp[NGROUPS*NSPECIES];

        
        
        // Initialize fluid velocities and Lorentz factor
        double n4D[4];
        double u4U[4];
        double v4U[4];
        double fidu_w_lorentz = 0.0;
        for(int a = 0; a < 4; a++) {
            if(a == 0) {
                u4U[a] = prims->u0;
                n4D[a] = -metric->lapse;
                v4U[a] = 0.0;
            }
            else {
                u4U[a] = prims->u0 * prims->vU[a-1];
                n4D[a] = 0;
                v4U[a] = prims->vU[a-1];
            }
            fidu_w_lorentz += -u4U[a] * n4D[a];
        }

        // Initialize projection tensor vector
        ghl_radiation_metric_tensor proj4;
        for(int a = 0; a < 4; a++) {
            for(int b = 0; b < 4; b++) {
            proj4.UD[a][b] = (int)(a == b);
                for(int c = 0; c < 4; c++) {
                    proj4.UD[a][b] += u4U[a] * u4U[c] * adm_aux->g4DD[b][c];
                }
            }
        }

        m1_root_params closure_params;
        closure_params.metric   = metric;
        closure_params.adm_aux  = adm_aux;
        closure_params.prims    = prims;
        closure_params.E        = 0.0;
        closure_params.F4       = &F4;
        closure_params.P4       = &P4;
        

        // Step 1 -- compute the sources
        printf("Step 1: Compute sources\n");
        for (int ig = 0; ig < NGROUPS*NSPECIES; ++ig) {
            if(M1_source_method == M1_source_method_explicit){
                //  Step exp-1-1 Initialize radiation variables
                printf("Step exp-1-1: Initialize radiation variables\n");
                E = rE_p[ig];
                F4.D[1] = rFx_p[ig];
                F4.D[2] = rFy_p[ig];
                F4.D[3] = rFz_p[ig];
                F4.D[0] = -metric->lapse * (n4D[1] * F4.D[1] +
                                        n4D[2] * F4.D[2] + 
                                        n4D[3] * F4.D[3]); 
                closure_params.E   = E;
                closure_params.F4  = &F4;
                closure_params.P4  = &P4;
                closure_params.chi = chi[ig];
                ghl_radiation_rootSolve_closure(&closure_params);
                

                double J = rJ[ig];
                double Gamma = compute_Gamma(
                    fidu_w_lorentz, &v4U, J, E, rad_E_floor, rad_eps, &F4);
                ghl_radiation_flux_vector H4 = {0};
                H4.D[0] = rHt[ig];
                H4.D[1] = rHx[ig];
                H4.D[2] = rHy[ig];
                H4.D[3] = rHz[ig];

                ghl_radiation_con_source_vector S4 = {0};
                calc_rad_sources(
                    eta_1[ig], abs_1[ig], scat_1[ig],
                    u4U, J, &H4, &S4);
                DrE[ig] = calc_rE_source(metric, &S4);
                ghl_radiation_con_source_vector rF_source = {0};
                calc_rF_source(metric, adm_aux, &S4, &rF_source);
                DrFx[ig] = rF_source.D[1];
                DrFy[ig] = rF_source.D[2];
                DrFz[ig] = rF_source.D[3];
                DrN[ig]  = cdt*metric->lapse*(volform*eta_0[ig] - abs_0[ig]*rN[ig]/Gamma);
            }
            
                                       
            // This is not done in THC: should I do closure here?
            
            // Step 1-1 Compute GR sources 
            printf("Step 1-1: Compute GR sources\n");
            // * This is done in thc_M1_calc_grsources.cc not in thc_M1_calc_update.cc
            // thc directly uses grid variables rE, rF to update GE, GF
            E = rE[ig];
            F4.D[1] = rFx[ig];
            F4.D[2] = rFy[ig];
            F4.D[3] = rFz[ig];
            F4.D[0] = -metric->lapse * (n4D[1] * F4.D[1] +
                                        n4D[2] * F4.D[2] + 
                                        n4D[3] * F4.D[3]);
            P4.DD[1][1] = rPxx[ig];
            P4.DD[1][2] = rPxy[ig];
            P4.DD[1][3] = rPxz[ig];
            P4.DD[2][1] = rPxy[ig];
            P4.DD[2][2] = rPyy[ig];
            P4.DD[2][3] = rPyz[ig];
            P4.DD[3][1] = rPxz[ig];
            P4.DD[3][2] = rPyz[ig];
            P4.DD[3][3] = rPzz[ig];

            double GE_source;
            ghl_radiation_con_source_vector GF_source;
            GE_source = calc_GE_source(&metric,
                           &metric_derivs_x,
                           &metric_derivs_y,
                           &metric_derivs_z,
                           &P4,
                           &F4,
                           &curv);
            calc_GF_source(&metric,
                           &metric_derivs_x,
                           &metric_derivs_y,
                           &metric_derivs_z,
                           E,
                           &F4,
                           &P4,
                           &GF_source);
            printf("GE_source = %f\n", GE_source);
            printf("GF_source: D = {%f, %f, %f, %f}\n",
                   GF_source.D[0], GF_source.D[1], GF_source.D[2], GF_source.D[3]);
            rE_rhs[ig]  = GE_source;
            rFx_rhs[ig] = GF_source.D[1];
            rFy_rhs[ig] = GF_source.D[2];
            rFz_rhs[ig] = GF_source.D[3];

            // Step 1-3 Compute 
            


            // TODO: Add a flag for only implement explicit update
            printf("Step 1-3: Advect radiation\n");
            //
            //
            // Here we boost to the fluid frame, compute fluid matter
            // interaction, and boost back. These values are used as
            // initial guess for the implicit solve.
            //
            // Advect radiation
            E_star = rE_p[ig] + cdt*rE_rhs[ig];
            
            F4_star.D[1] = rFx_p[ig] + cdt*rFx_rhs[ig];
            F4_star.D[2] = rFy_p[ig] + cdt*rFy_rhs[ig];
            F4_star.D[3] = rFz_p[ig] + cdt*rFz_rhs[ig];
            // double rad_eps = 1.0e-15; // default value in THC
            // double rad_E_floor = 1.0e-15;
            // apply_floor(adm_aux, &E, &F4, rad_E_floor, rad_eps); // TODO: what does this step do?

            apply_floor(adm_aux, &E_star, &F4_star, rad_E_floor, rad_eps);
            double N_star = fmax(rN_p[ig] + cdt*rN_rhs[ig], rad_N_floor);

            //
            // Compute quantities in the fluid frame
            
            closure_params.E        = E_star;
            closure_params.F4       = &F4_star;
            closure_params.P4       = &P4;
            ghl_radiation_rootSolve_closure(&closure_params);
            // Need to calculate if P4 is inf

            
            assemble_rT_lab_frame(n4D, E_star, &F4_star, &P4, &rT4DD);

            double const Jstar = calc_J_from_rT(u4U, &proj4, &rT4DD);

            ghl_radiation_flux_vector H4_star = {0};
            calc_H4D_from_rT(u4U, &proj4, &rT4DD, &H4_star);

            printf("Step 1-4: Estimate interaction with matter\n");
            //
            // Estimate interaction with matter
            double const dtau = cdt/fidu_w_lorentz;
            double J_new = (Jstar + dtau*eta_1[ig]*volform)/(1 + dtau*abs_1[ig]);

            // Only three components of H^a are independent H^0 is found by
            // requiring H^a u_a = 0
            double const khat = (abs_1[ig] + scat_1[ig]);

            ghl_radiation_flux_vector H4_new = {0};
            for (int a = 1; a < 4; ++a) {
                H4_new.D[a] = H4_star.D[a]/(1 + dtau*khat);
            }
            H4_new.D[0] = 0.0;
            for (int a = 1; a < 4; ++a) {
                H4_new.D[0] -= H4_new.D[a]*(u4U[a]/u4U[0]);
            }

            printf("Step 1-5: Update Tmunu\n");
            //
            // Update Tmunu
            double const H2 = ghl_compute_vec2_from_vec4D(adm_aux->g4UU, H4_new.D);
            double const dthick = 3.*(1. - chi[ig])/2.;
            double const dthin = 1. - dthick;

            for(int a = 0; a < 4; ++a) {
                for(int b = 0; b < 4; ++b) {
                    rT4DD.T4[a][b] = J_new*u4U[a]*u4U[b] 
                                    + H4_new.D[a]*u4U[b] + H4_new.D[b]*u4U[a] +
                        dthin*J_new*(H4_new.D[a]*H4_new.D[b]*(H2 > 0 ? 1/H2 : 0)) +
                        dthick*J_new*(adm_aux->g4DD[a][b] + u4U[a]*u4U[b])/3;
                }
            }

            printf("Step 1-6: Source update\n");
            //
            // Boost back to the lab frame
            // Initialize E_new, F4_new for source update
            E_new = calc_J_from_rT(u4U, &proj4, &rT4DD);
            calc_H4D_from_rT(u4U, &proj4, &rT4DD, &F4_star);
            apply_floor(adm_aux, &E_new, &F4_new, rad_E_floor, rad_eps);
            //
            // Compute interaction with matter
            ghl_m1_powell_params p     = {0}; 
            ghl_radiation_con_source_vector rF_source_0 = {0};

            // int source_err = ghl_source_update(thc_params,&p, &E_new, &F4_new);
            int source_err = ghl_source_update(
                &thc_params, chi[ig], eta_0[ig], abs_0[ig], scat_1[ig], cdt,
                metric, adm_aux, prims,
                E_star, &F4_star,
                closure, gsl_solver_1d, gsl_solver_nd,
                &E_new, &F4_new);
            int souce_err = 0;
            apply_floor(adm_aux, &E_new, &F4_new, rad_E_floor, rad_eps);


            //
            // Update closure
            ghl_radiation_rootSolve_closure(&closure_params);
            // ghl_radiation_apply_closure(metric, adm_aux, prims, E_new, &F4_new, chi[ig], &P4);

            printf("Step 1-7: Compute radiation changes\n");
            //
            // Compute new radiation energy density in the fluid frame
            ghl_stress_energy T4;
            assemble_rT_lab_frame(n4D, E_new, &F4_new, &P4, &T4);
            J_new = calc_J_from_rT(u4U, &proj4, &T4);

            
            //
            // Compute changes in radiation energy and momentum
            DrE[ig]  = E_new - E_star;
            DrFx[ig] = F4_new.D[1] - F4_star.D[1];
            DrFy[ig] = F4_new.D[2] - F4_star.D[2];
            DrFz[ig] = F4_new.D[3] - F4_star.D[3];

            //
            // Compute updated Gamma
            double const Gamma = compute_Gamma(fidu_w_lorentz, &v4U, J_new, E_new, rad_E_floor, rad_eps, &F4_new);

            //
            // N^k+1 = N^* + dt ( eta - abs N^k+1 )
            if (source_therm_limit < 0 || cdt*abs_0[ig] < source_therm_limit) {
                DrN[ig] = (N_star + cdt*metric->lapse*volform*eta_0[ig])/
                            (1 + cdt*metric->lapse*abs_0[ig]/Gamma) - N_star;
            }
            //
            // The neutrino number density is updated assuming the neutrino
            // average energies are those of the equilibrium
            // TODO: need to calculate this as calc_Opacity in thc
            else {
                DrN[ig] = (nueave[ig] > 0 ? Gamma*J_new/nueave[ig] - N_star : 0.0);
            }

            //
            // Fluid lepton sources
            DDxp[ig] = -mb*(DrN[ig]*(ig == 0) - DrN[ig]*(ig == 1));
        } // for (int ig ...)

        //
        // Step 2 -- limit the sources
        printf("Step 2: Limit sources\n");
        double theta = 1.0;
        

        if (source_limiter >= 0) {
            theta = 1.0;
            double DTau_sum = 0.0;
            double DDxp_sum = 0.0;
            for (int ig = 0; ig < NGROUPS*NSPECIES; ++ig) {
                double E_star = rE_p[ig] + cdt*rE_rhs[ig];
                if (DrE[ig] < 0) {
                    theta = fmin(-source_limiter*fmax(E_star, 0.0)/DrE[ig], theta);
                }
                DTau_sum -= DrE[ig];

                double N_star = rN_p[ig] + cdt*rN_rhs[ig];
                if (DrN[ig] < 0) {
                    theta = fmin(-source_limiter*fmax(N_star, 0.0)/DrN[ig], theta);
                }
                DDxp_sum += DDxp[ig];
            }
            double const DYe = DDxp_sum/dens;
            if (DTau_sum < 0) {
                theta = fmin(-source_limiter*fmax(tau, 0.0)/DTau_sum, theta);
            }
            if (DYe > 0) {
                theta = fmin(source_limiter*fmax(source_Ye_max - Y_e, 0.0)/DYe, theta);
            }
            else if (DYe < 0) {
                theta = fmin(source_limiter*fmin(source_Ye_min - Y_e, 0.0)/DYe, theta);
            }
            theta = fmax(0.0, theta);
        }

        printf("Step 3: Update radiation fields\n");
        //
        // Step 3 -- update fields
        for (int ig = 0; ig < NGROUPS*NSPECIES; ++ig) {
            //
            // Update radiation quantities
            double E     =  rE_p[ig] + cdt*rE_rhs[ig]  + theta*DrE[ig];
            F4.D[1]      = rFx_p[ig] + cdt*rFx_rhs[ig] + theta*DrFx[ig];
            F4.D[2]      = rFy_p[ig] + cdt*rFy_rhs[ig] + theta*DrFy[ig];
            F4.D[3]      = rFz_p[ig] + cdt*rFz_rhs[ig] + theta*DrFz[ig];
            apply_floor(adm_aux, &E, &F4, rad_E_floor, rad_eps);

            double N =  rN_p[ig] + cdt*rN_rhs[ig]  + theta*DrN[ig];
            N = fmax(N, rad_N_floor);

            // TODO: comment out this section and work on it later
            // Compute back reaction on the fluid
            // NOTE: fluid backreaction is only needed at the last substep
            // if (backreact && 0 == *TimeIntegratorStage) {
            //     assert (NGROUPS == 1);
            //     assert (NSPECIES == 3);

            //     sconx  -= theta*DrFx[ig];
            //     scony  -= theta*DrFy[ig];
            //     sconz  -= theta*DrFz[ig];
            //     tau    -= theta*DrE[ig];
            //     densxp += theta*DDxp[ig];
            //     densxn -= theta*DDxp[ig];

            //     netabs  += theta*DDxp[ig];
            //     netheat -= theta*DrE[ig];
            // }

            //
            // Save updated results into grid functions
            rE[ig]  = E;
            rFx[ig] = F4.D[1];
            rFy[ig] = F4.D[2];
            rFz[ig] = F4.D[3];
            rN[ig]  = N;
            printf("Results: \n");
            printf("E: %f  F4: (%f, %f, %f, %f)  N: %f\n",
                   rE[ig], rFx[ig], rFy[ig], rFz[ig], rN[ig]);
        }
        gsl_root_fsolver_free(gsl_solver_1d);
        gsl_multiroot_fdfsolver_free(gsl_solver_nd);
        // End parallel 

        // Restore GSL error handler
        gsl_set_error_handler(gsl_err);
    } // for (int timeIntegrator ...)
} // ghl_M1_update