#include "ghl.h"
#include "ghl_m1.h"
#include "ghl_radiation.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

// This file translates from thc_M1_calc_update.cc
#define NGROUPS 1
#define NSPECIES 1
int ghl_m1_update(
    double timestep,
    const ghl_metric_quantities *metric,
    const ghl_extrinsic_curvature *curv,
    const ghl_ADM_aux_quantities *adm_aux,
    const ghl_primitive_quantities *prims,
    const ghl_m1_thc_params *thc_params,
    const double abs_0[NGROUPS*NSPECIES],
    const double abs_1[NGROUPS*NSPECIES],
    const double eta_0[NGROUPS*NSPECIES],
    const double eta_1[NGROUPS*NSPECIES],
    const double scat_0[NGROUPS*NSPECIES],
    const double scat_1[NGROUPS*NSPECIES],
    ghl_m1_conservative_quantities *rU,
    ghl_m1_conservative_quantities *rU_p,
    ghl_m1_conservative_quantities *rhs_flux,
    ghl_m1_conservative_quantities *rhs_gr
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
    
    ///////////////////////////////////

    
    // Public grid variables (in interface.ccl)
    // Evolved quantities in current time step
    // thc_M1_setsym.cc        :: Setup
    double  rE[NGROUPS*NSPECIES]   = {0};
    double rFx[NGROUPS*NSPECIES]   = {0};
    double rFy[NGROUPS*NSPECIES]   = {0};
    double rFz[NGROUPS*NSPECIES]   = {0};
    double  rN[NGROUPS*NSPECIES]   = {0};

    // Evolved quantities in previous time step
    // thc_M1_copy_levels.cc:    :: std::memcpy(rE_p, rE, siz); copy rE to rE_p
    // thc_M1_initial_copy.cc:   :: std::memcpy(rE, rE_p, siz); copy rE_p to rE
    
    // Read
    double  rE_p[NGROUPS*NSPECIES] = {0};
    double rFx_p[NGROUPS*NSPECIES] = {0};
    double rFy_p[NGROUPS*NSPECIES] = {0};
    double rFz_p[NGROUPS*NSPECIES] = {0};
    double  rN_p[NGROUPS*NSPECIES] = {0};

    // thc_M1_init_rhs.cc       :: First defined here with size = nspecies * sizeof(CCTK_REAL)
    // thc_M1_calc_flux.cc      :: append rhs pointer to rhs grid variable in  (rhs(ig,0) = &rE_rhs[i4D])
    //                             add flux term to rhs
    // thc_M1_calc_grsources.cc :: add source term to rE_rhs, rF_rhs
    // thc_M1_calc_update.cc    :: Read rE_rhs, rF_rhs to update Estar, Fstar (or directly apply for explicit)
    double rE_rhs[NGROUPS*NSPECIES]   = {0};
    double rFx_rhs[NGROUPS*NSPECIES]  = {0};
    double rFy_rhs[NGROUPS*NSPECIES]  = {0};
    double rFz_rhs[NGROUPS*NSPECIES]  = {0};
    double rN_rhs[NGROUPS*NSPECIES]   = {0};
    
    // From function arguments store to local arrays
    // these variables will not be changed until update is done
    for (int ig = 0; ig < NGROUPS*NSPECIES; ++ig) {
        rE[ig]  = rU_p->rE[ig];
        rFx[ig] = rU_p->rF4[ig].D[1];
        rFy[ig] = rU_p->rF4[ig].D[2];
        rFz[ig] = rU_p->rF4[ig].D[3];
        rN[ig]  = rU_p->rN[ig];

        rE_p[ig]  = rU_p->rE[ig];
        rFx_p[ig] = rU_p->rF4[ig].D[1];
        rFy_p[ig] = rU_p->rF4[ig].D[2];
        rFz_p[ig] = rU_p->rF4[ig].D[3];
        rN_p[ig]  = rU_p->rN[ig];

        // RHS
        // Initialize the RHS to zero
        rE_rhs[ig]  = 0.0;
        rFx_rhs[ig] = 0.0;
        rFy_rhs[ig] = 0.0;
        rFz_rhs[ig] = 0.0;
        rN_rhs[ig]  = 0.0;

        // GR sources
        rE_rhs[ig]  += rhs_gr->rE[ig];
        rFx_rhs[ig] += rhs_gr->rF4[ig].D[1];
        rFy_rhs[ig] += rhs_gr->rF4[ig].D[2];
        rFz_rhs[ig] += rhs_gr->rF4[ig].D[3];
        rN_rhs[ig]  += rhs_gr->rN[ig];
        // printf("rE_rhs[%d] after GR source = %f\n", ig, rE_rhs[ig]);

        // fluxes
        rE_rhs[ig]  += rhs_flux->rE[ig];
        rFx_rhs[ig] += rhs_flux->rF4[ig].D[1];
        rFy_rhs[ig] += rhs_flux->rF4[ig].D[2];
        rFz_rhs[ig] += rhs_flux->rF4[ig].D[3];
        rN_rhs[ig]  += rhs_flux->rN[ig];
        // printf("rE_rhs[%d] after flux source = %f\n", ig, rE_rhs[ig]);
    }

    //Read
    double nueave[NSPECIES*NGROUPS]   = {0};
    double chi[NGROUPS*NSPECIES]      = {0};
        
    // opacities

    //Read
    // These are probably from calc_opacity, do they exist in leakage already?
    // double abs_0[NGROUPS*NSPECIES]    = {0};
    // double abs_1[NGROUPS*NSPECIES]    = {0};
    // double eta_0[NGROUPS*NSPECIES]    = {0};
    // double eta_1[NGROUPS*NSPECIES]    = {0};
    // double scat_1[NGROUPS*NSPECIES]   = {0};

    // THC flag
    int M1_source_method = 1; // 0: explicit, 1: implicit, 2: boost
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

    // In WhiskyTHC, thcextra/Weakrates/atmois_mass
    // AtomicMassImpl = normfact * cgs2cactusMass * mass_fact * mev_to_erg / (clight*clight)
    double const mb = 1.0; //AverageBaryonMass(); //TODO: what exactly is this function?
    
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

        // write
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
        // printf("Step 1: Compute sources\n");
        for (int ig = 0; ig < NGROUPS*NSPECIES; ++ig) {
            // Let's hope we never use this
            if(M1_source_method == M1_source_method_explicit) {
                //  Step exp-1-1 Initialize radiation variables
                // printf("Step exp-1-1: Initialize radiation variables\n");

                // In THC this closure step is not needed because J and H are read from grid variables
                // but we don't have those so we compute them here
                // E, F are needed for Gamma;  J, H are needed for S4
                double E = rE_p[ig];
                ghl_radiation_flux_vector F4 = {0};
                ghl_radiation_pressure_tensor P4 = {0};
                F4.D[1] = rFx_p[ig];
                F4.D[2] = rFy_p[ig];
                F4.D[3] = rFz_p[ig];
                // F4.D[0] = -metric->lapse * (n4D[1] * F4.D[1] +
                //                         n4D[2] * F4.D[2] + 
                //                         n4D[3] * F4.D[3]); 
                F4.D[0] = metric->betaU[1] * F4.D[1] + metric->betaU[2] * F4.D[2] + metric->betaU[3] * F4.D[3];
                closure_params.E   = E;
                closure_params.F4  = &F4;
                closure_params.P4  = &P4;
                closure_params.chi = chi[ig];
                ghl_radiation_rootSolve_closure(&closure_params);
                
                ghl_stress_energy T4_exp = {0};
                assemble_rT_lab_frame(n4D, E, &F4, &P4, &T4_exp);
                double J = calc_J_from_rT(u4U, &proj4, &T4_exp);
                ghl_radiation_flux_vector H4 = {0};
                calc_H4D_from_rT(u4U, &proj4, &T4_exp, &H4);
                // printf("J = %f, H4 = {%f, %f, %f, %f}\n", J, H4.D[0], H4.D[1], H4.D[2], H4.D[3]);
                double Gamma = compute_Gamma(
                    fidu_w_lorentz, &v4U, J, E, rad_E_floor, rad_eps, &F4);
                

                ghl_radiation_con_source_vector S4 = {0};
                // printf("eta_1 = %f, abs_1 = %f, scat_1 = %f\n", eta_1[ig], abs_1[ig], scat_1[ig]);
                calc_rad_sources(
                    eta_1[ig], abs_1[ig], scat_1[ig],
                    u4U, J, &H4, &S4);
                // printf("S4={%f, %f, %f, %f}\n", S4.U[0], S4.U[1], S4.U[2], S4.U[3]);
                DrE[ig] = calc_rE_source(metric, &S4);

                ghl_radiation_con_source_vector rF_source = {0};
                calc_rF_source(metric, adm_aux, &S4, &rF_source);
                DrFx[ig] = cdt * rF_source.D[1];
                DrFy[ig] = cdt * rF_source.D[2];
                DrFz[ig] = cdt * rF_source.D[3];
                DrN[ig]  = cdt*metric->lapse*(volform*eta_0[ig] - abs_0[ig]*rN[ig]/Gamma);
            } else { // not explicit
                // Step 1-1 Compute GR sources 
                // printf("Step 1-1: Compute GR sources\n");
                // * This is done in thc_M1_calc_grsources.cc not in thc_M1_calc_update.cc
                // thc directly uses grid variables rE, rF to update GE, GF
                rE_rhs[ig]  = rhs_gr->rE[ig];
                rFx_rhs[ig] = rhs_gr->rF4[ig].D[1];
                rFy_rhs[ig] = rhs_gr->rF4[ig].D[2];
                rFz_rhs[ig] = rhs_gr->rF4[ig].D[3];
                rN_rhs[ig]  = rhs_gr->rN[ig];
                // printf("grsources: rhs_gr->rE[%d] = %f, rhs_gr->rF4[%d].D[1] = %f, rhs_gr->rF4[%d].D[2] = %f, rhs_gr->rF4[%d].D[3] = %f, rhs_gr->rN[%d] = %f\n", ig, rhs_gr->rE[ig], ig, rhs_gr->rF4[ig].D[1], ig, rhs_gr->rF4[ig].D[2], ig, rhs_gr->rF4[ig].D[3], ig, rhs_gr->rN[ig]);

                // Step 1-2 Compute flux sources (not done in THC)
                // printf("Step 1-2: Compute flux sources\n");
                rE_rhs[ig]  += rhs_flux->rE[ig];
                rFx_rhs[ig] += rhs_flux->rF4[ig].D[1];
                rFy_rhs[ig] += rhs_flux->rF4[ig].D[2];
                rFz_rhs[ig] += rhs_flux->rF4[ig].D[3];
                rN_rhs[ig]  += rhs_flux->rN[ig];
                // printf("flux sources: rhs_flux->rE[%d] = %f, rhs_flux->rF4[%d].D[1] = %f, rhs_flux->rF4[%d].D[2] = %f, rhs_flux->rF4[%d].D[3] = %f, rhs_flux->rN[%d] = %f\n", ig, rhs_flux->rE[ig], ig, rhs_flux->rF4[ig].D[1], ig, rhs_flux->rF4[ig].D[2], ig, rhs_flux->rF4[ig].D[3], ig, rhs_flux->rN[ig]);
                

                // Step 1-3 Advect radiation (not done in THC)
                // printf("Step 1-3: Advect radiation\n");
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

                apply_floor(adm_aux, &E_star, &F4_star, rad_E_floor, rad_eps); 
                // printf("E_star[%d] = %f, F4_star[%d] = {%f, %f, %f}\n", ig, E_star, ig, F4_star.D[1], F4_star.D[2], F4_star.D[3]);  
                
                double N_star = fmax(rN_p[ig] + cdt*rN_rhs[ig], rad_N_floor);

                //
                // Compute quantities in the fluid frame
                
                closure_params.E        = E_star;
                closure_params.F4       = &F4_star;
                closure_params.P4       = &P4;
                ghl_radiation_rootSolve_closure(&closure_params);
                assemble_rT_lab_frame(n4D, E_star, &F4_star, &P4, &rT4DD);
                // Need to calculate if P4 is inf

                double const Jstar = calc_J_from_rT(u4U, &proj4, &rT4DD);
                ghl_radiation_flux_vector H4_star = {0};
                calc_H4D_from_rT(u4U, &proj4, &rT4DD, &H4_star);
                // printf("Jstar = %f, H4_star = {%f, %f, %f, %f}\n", Jstar, H4_star.D[0], H4_star.D[1], H4_star.D[2], H4_star.D[3]);

                // printf("Step 1-4: Estimate interaction with matter\n");
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

                // printf("Step 1-5: Update Tmunu\n");
                //
                // Update Tmunu
                double const H2 = ghl_compute_vec2_from_vec4D(adm_aux->g4UU, H4_new.D);
                // if (M1_SRC_METHOD == "boost") {
                //     CCTK_REAL const xi = sqrt(H2)*(Jnew > rad_E_floor ? 1/Jnew : 0);
                //     chi[i4D] = closure_fun(xi);
                // }
                // this is really only important in the thick limit, so take chi = 1/3
                chi[ig] = 1.0/3; 

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
                // printf("updated Tmunu: rT4DD.T4 = \n{%f, %f, %f, %f},\n{%f, %f, %f, %f},\n{%f, %f, %f, %f},\n{%f, %f, %f, %f}\n",
                //     rT4DD.T4[0][0], rT4DD.T4[0][1], rT4DD.T4[0][2], rT4DD.T4[0][3],
                //     rT4DD.T4[1][0], rT4DD.T4[1][1], rT4DD.T4[1][2], rT4DD.T4[1][3],
                //     rT4DD.T4[2][0], rT4DD.T4[2][1], rT4DD.T4[2][2], rT4DD.T4[2][3],
                //     rT4DD.T4[3][0], rT4DD.T4[3][1], rT4DD.T4[3][2], rT4DD.T4[3][3]);

                // printf("Step 1-6: Source update\n");
                //
                // Boost back to the lab frame
                // Initialize E_new, F4_new for source update
                double n4U[4] = {metric->lapseinv, 
                        -metric->betaU[0] * metric->lapseinv,
                        -metric->betaU[1] * metric->lapseinv,
                        -metric->betaU[2] * metric->lapseinv };
                E_new = calc_J_from_rT(n4U, &proj4, &rT4DD);
                // calc_H4D_from_rT(u4U, &proj4, &rT4DD, &H4_star);
                calc_H4D_from_rT(n4U, &proj4, &rT4DD, &F4_new);
                apply_floor(adm_aux, &E_new, &F4_new, rad_E_floor, rad_eps); 
                // printf("Initial guess E_new[%d] = %f, F4_new[%d] = {%f, %f, %f}\n", ig, E_new, ig, F4_new.D[1], F4_new.D[2], F4_new.D[3]);  

                //
                // Compute interaction with matter
                // TODO: What is going on here:
                //ghl_m1_powell_params p     = {0}; 
                //ghl_radiation_con_source_vector rF_source_0 = {0};

                // int source_err = ghl_source_update(thc_params,&p, &E_new, &F4_new);
                int source_err = ghl_source_update(
                    &thc_params, chi[ig], eta_1[ig], abs_1[ig], scat_1[ig], cdt,
                    metric, adm_aux, prims,
                    E_star, &F4_star,
                    closure, gsl_solver_1d, gsl_solver_nd,
                    &E_new, &F4_new);
                int souce_err = 0;
                apply_floor(adm_aux, &E_new, &F4_new, rad_E_floor, rad_eps); 


                //
                // Update closure
                closure_params.E  = E_new;
                closure_params.F4 = &F4_new;
                closure_params.P4 = &P4;
                ghl_radiation_rootSolve_closure(&closure_params);
                // ghl_radiation_apply_closure(metric, adm_aux, prims, E_new, &F4_new, chi[ig], &P4);

                // printf("Step 1-7: Compute radiation changes\n");
                //
                // Compute new radiation energy density in the fluid frame
                ghl_stress_energy T4;
                assemble_rT_lab_frame(n4D, E_new, &F4_new, &P4, &T4);
                // printf("Tmunu: rT4DD.T4 = \n{%f, %f, %f, %f},\n{%f, %f, %f, %f},\n{%f, %f, %f, %f},\n{%f, %f, %f, %f}\n",
                //     T4.T4[0][0], T4.T4[0][1], T4.T4[0][2], T4.T4[0][3],
                //     T4.T4[1][0], T4.T4[1][1], T4.T4[1][2], T4.T4[1][3],
                //     T4.T4[2][0], T4.T4[2][1], T4.T4[2][2], T4.T4[2][3],
                //     T4.T4[3][0], T4.T4[3][1], T4.T4[3][2], T4.T4[3][3]);
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
            } // end if not explicit

            //
            // Fluid lepton sources
            DDxp[ig] = -mb*(DrN[ig]*(ig == 0) - DrN[ig]*(ig == 1));
        } // for (int ig ...)
        // printf("E: %f, F4: {%f, %f, %f, %f}, N: %f\n", E_new, F4_new.D[0], F4_new.D[1], F4_new.D[2], F4_new.D[3], 0.0); // N not included

        //
        // Step 2 -- limit the sources
        // printf("Step 2: Limit sources\n");
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
        // printf("theta = %f\n", theta);

        // printf("Step 3: Update radiation fields\n");
        // //
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
        } // for (int ig ...)
        gsl_root_fsolver_free(gsl_solver_1d);
        gsl_multiroot_fdfsolver_free(gsl_solver_nd);
        // End parallel 

        // Restore GSL error handler
        gsl_set_error_handler(gsl_err);
    } // for (int timeIntegrator ...)
    for (int ig = 0; ig < NGROUPS*NSPECIES; ++ig) {
        //Store previous time step
        rU_p->rE[ig] = rU->rE[ig];
        rU_p->rF4[ig].D[1] = rU->rF4[ig].D[1];
        rU_p->rF4[ig].D[2] = rU->rF4[ig].D[2];
        rU_p->rF4[ig].D[3] = rU->rF4[ig].D[3];
        rU_p->rN[ig] = rU->rN[ig];

        //update time step
        rU->rE[ig] = rE[ig];
        rU->rF4[ig].D[1] = rFx[ig];
        rU->rF4[ig].D[2] = rFy[ig];
        rU->rF4[ig].D[3] = rFz[ig];
        rU->rN[ig] = rN[ig];
    }
  return 0;
} // ghl_M1_update