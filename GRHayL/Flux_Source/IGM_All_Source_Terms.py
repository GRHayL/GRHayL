# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import os,sys
import GRMHD_equations_new_version as GRMHD    # NRPy+: Generate general relativistic magnetohydrodynamics equations

nrpy_dir_path = os.path.join("nrpy/")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)

# Step P1: Import needed NRPy+ core modules:
from outputC import outputC, outCfunction # NRPy+: Core C code output module
import NRPy_param_funcs as par        # NRPy+: Parameter interface
import reference_metric as rfm        # NRPy+: Reference metric support
import sympy as sp                # SymPy: The Python computer algebra package upon which NRPy+ depends

thismodule = __name__

par.initialize_param(par.glb_param(type="bool", module=thismodule, 
       parname="using_Valencia_velocity", defaultval=True))

#Step 0: Set the spatial dimension parameter to 3.
par.set_parval_from_str("grid::DIM", 3)
DIM = par.parval_from_str("grid::DIM")

CoordSystem = "Cartesian"

# Set coordinate system to dst_basis
par.set_parval_from_str("reference_metric::CoordSystem", CoordSystem)
rfm.reference_metric()

def Cfunction__GRMHD_SourceTerms(Ccodesdir, includes=None, formalism="ADM", outCparams = "outCverbose=False,CSE_sorting=canonical,CSE_enable=True"):
    # Generate SymPy symbolic expressions
    GRMHD.set_up_base_vars(formalism=formalism)

    GRMHD.compute_vU_from_u4U__no_speed_limit(GRMHD.u4U)   

    GRMHD.compute_sqrtgammaDET(GRMHD.gammaDD)
    GRMHD.compute_smallb4U(GRMHD.gammaDD,GRMHD.betaU,GRMHD.alpha, GRMHD.u4U, GRMHD.BU, GRMHD.sqrt4pi)
    GRMHD.compute_smallbsquared(GRMHD.gammaDD,GRMHD.betaU,GRMHD.alpha, GRMHD.smallb4U)

    # First compute stress-energy tensor T4UU and T4UD:
    GRMHD.compute_T4UU(GRMHD.gammaDD,GRMHD.betaU,GRMHD.alpha, GRMHD.rho_b,GRMHD.P,GRMHD.h,GRMHD.u4U, GRMHD.smallb4U, GRMHD.smallbsquared)
    GRMHD.compute_T4UD(GRMHD.gammaDD,GRMHD.betaU,GRMHD.alpha, GRMHD.T4UU)

    # Compute conservative variables in terms of primitive variables
    GRMHD.compute_rho_star(GRMHD.alpha, GRMHD.sqrtgammaDET, GRMHD.rho_b,GRMHD.u4U)
    GRMHD.compute_tau_tilde(GRMHD.alpha, GRMHD.sqrtgammaDET, GRMHD.T4UU,GRMHD.rho_star)
    GRMHD.compute_S_tildeD(GRMHD.alpha, GRMHD.sqrtgammaDET, GRMHD.T4UD)

    # Next compute fluxes of conservative variables
    GRMHD.compute_rho_star_fluxU(GRMHD.VU,GRMHD.rho_star)
    GRMHD.compute_tau_tilde_fluxU(GRMHD.alpha, GRMHD.sqrtgammaDET, GRMHD.VU, GRMHD.T4UU, GRMHD.rho_star)
    GRMHD.compute_S_tilde_fluxUD (GRMHD.alpha, GRMHD.sqrtgammaDET, GRMHD.T4UD)

    # Then declare derivatives & compute g4DD_zerotimederiv_dD
    GRMHD.compute_g4DD_zerotimederiv_dD(GRMHD.gammaDD,GRMHD.betaU,GRMHD.alpha, GRMHD.gammaDD_dD,GRMHD.betaU_dD,GRMHD.alpha_dD)

    # Then compute source terms on tau_tilde and S_tilde equations
    GRMHD.compute_tau_tilde_source_term(GRMHD.KDD,GRMHD.betaU,GRMHD.alpha, GRMHD.sqrtgammaDET, GRMHD.alpha_dD, GRMHD.T4UU)
    GRMHD.compute_S_tilde_source_termD(GRMHD.alpha,GRMHD.sqrtgammaDET,GRMHD.g4DD_zerotimederiv_dD, GRMHD.T4UU)
    
    # tau_tilde_source_term_free_symbols = GRMHD.tau_tilde_source_term.free_symbols
    S_tilde_source_termD0_free_symbols = GRMHD.S_tilde_source_termD[0].free_symbols
    S_tilde_source_termD1_free_symbols = GRMHD.S_tilde_source_termD[1].free_symbols
    S_tilde_source_termD2_free_symbols = GRMHD.S_tilde_source_termD[2].free_symbols

    # all_free_sysmbols = tau_tilde_source_term_free_symbols.union(\
    all_free_sysmbols = S_tilde_source_termD0_free_symbols.union(\
                                 S_tilde_source_termD1_free_symbols,
                                 S_tilde_source_termD2_free_symbols)

    prims_NRPy = ["u4U0", "u4U1", "u4U2", "u4U3", "BU0", "BU1", "BU2", "P", "rhob"]
    prims_GRHayL = ["u0", "vx*u4U0", "vy*u4U0", "vz*u4U0", "Bx", "By", "Bz", "press", "rho"]

    prestring = r"""
double h, cs2;

eos->compute_h_and_cs2(eos, prims, &h, &cs2);
"""

    for i in range(len(prims_NRPy)):
        prestring += "const double "+prims_NRPy[i]+" = prims->"+prims_GRHayL[i]+";\n"

    prestring += "const double "+str(GRMHD.alpha)+" = metric->lapse;\n"

    checker = []


    if formalism=="BSSN":
        # BSSN quantites
        prestring += "const double "+str(GRMHD.Bq.trK)+" = metric_quantities->"+str(GRMHD.Bq.trK)+";\n"
        prestring += "const double "+str(GRMHD.Bq.cf)+" = metric_quantities->"+str(GRMHD.Bq.cf)+";\n"

        for i in range(3):
            vetU_var = GRMHD.Bq.vetU[i]
            prestring += "const double "+str(vetU_var)+" = metric_quantities->"+str(vetU_var)+";\n"

        for i in range(3):
            for j in range(3):
                aDD_var = GRMHD.Bq.aDD[i][j]
                if aDD_var in checker: 
                    continue
                prestring += "const double "+str(aDD_var)+" = metric_quantities->"+str(aDD_var)+";\n"
                checker.append(aDD_var)

        for i in range(3):
            for j in range(3):
                hDD_var = GRMHD.Bq.hDD[i][j]
                if hDD_var in checker: 
                    continue
                prestring += "const double "+str(hDD_var)+" = metric_quantities->"+str(hDD_var)+";\n"
                checker.append(hDD_var)

        for var in all_free_sysmbols:
            if "_dD" in str(var):
                prestring += "const double "+str(var)+" = metric_quantities_derivatives->"+str(var)+";\n"

    else:
        #ADM quantities
    #     for i in range(3):
    #             betaU_var = GRMHD.betaU[i]
    #             prestring += "const double "+str(betaU_var)+" = metric_quantities->"+str(betaU_var)+";\n"

    #     for i in range(3):
    #         for j in range(3):
    #             KDD_var = GRMHD.KDD[i][j]
    #             if KDD_var in checker: 
    #                 continue
    #             prestring += "const double "+str(KDD_var)+" = metric_quantities->"+str(KDD_var)+";\n"
    #             checker.append(KDD_var)

    #     for i in range(3):
    #         for j in range(3):
    #             gammaDD_var = GRMHD.gammaDD[i][j]
    #             if gammaDD_var in checker: 
    #                 continue                
    #             prestring += "const double "+str(gammaDD_var)+" = metric_quantities->"+str(gammaDD_var)+";\n"
    #             checker.append(gammaDD_var)

    #     for var in all_free_sysmbols:
    #         if "_dD" in str(var):
    #             prestring += "const double "+str(var)+" = metric_quantities_derivatives->"+str(var)+";\n"


        prestring += "const double betaU0 = metric->betax;\n"
        prestring += "const double betaU1 = metric->betay;\n"
        prestring += "const double betaU2 = metric->betaz;\n"


        prestring += "const double gammaDD00 = metric->adm_gxx;\n"
        prestring += "const double gammaDD01 = metric->adm_gxy;\n"
        prestring += "const double gammaDD02 = metric->adm_gxz;\n"

        prestring += "const double gammaDD11 = metric->adm_gyy;\n"
        prestring += "const double gammaDD12 = metric->adm_gyz;\n"

        prestring += "const double gammaDD22 = metric->adm_gzz;\n"
    
        vars_to_write = ["cons->S_x", "cons->S_y", "cons->S_z", "cons->tau"]
        
        vars_rhs = [GRMHD.S_tilde_source_termD[0], 
                    GRMHD.S_tilde_source_termD[1], 
                    GRMHD.S_tilde_source_termD[2], 
                    GRMHD.tau_tilde_source_term_2_3_A, 
                    GRMHD.tau_tilde_source_term_2_3_B, 
                    GRMHD.tau_tilde_source_term_2_3_C,
                    GRMHD.tau_tilde_source_term_1_3,
                    ]

        c_type = "void"

        params  = "const primitive_quantities *restrict prims, "
        params  += "struct eos_parameters const *restrict eos, "
        params  += "const metric_quantities *restrict metric, "
        params  += "const metric_derivatives *restrict metric_derivs, "
        params  += "conservative_quantities *restrict cons"

        for i in range(3):
            desc = f"Add source term for {i}-component of Stilde and tau_tilde"
            name = f"calculate_source_terms_dirn{i}"

            loopstring = prestring
            loopstring += f"const double alpha_dD{i} = metric_derivs->lapse[{i}];\n"

            loopstring += f"const double betaU_dD0{i} = metric_derivs->betax[{i}];\n"
            loopstring += f"const double betaU_dD1{i} = metric_derivs->betay[{i}];\n"
            loopstring += f"const double betaU_dD2{i} = metric_derivs->betaz[{i}];\n"

            loopstring += f"const double gammaDD_dD00{i} = metric_derivs->adm_gxx[{i}];\n"
            loopstring += f"const double gammaDD_dD01{i} = metric_derivs->adm_gxy[{i}];\n"
            loopstring += f"const double gammaDD_dD02{i} = metric_derivs->adm_gxz[{i}];\n"

            loopstring += f"const double gammaDD_dD11{i} = metric_derivs->adm_gyy[{i}];\n"
            loopstring += f"const double gammaDD_dD12{i} = metric_derivs->adm_gyz[{i}];\n"

            loopstring += f"const double gammaDD_dD22{i} = metric_derivs->adm_gzz[{i}];\n"


            body = outputC([vars_rhs[i], vars_rhs[i+3]], [vars_to_write[i], vars_to_write[-1]],
                                  params=outCparams, 
                       filename="returnstring", prestring=loopstring)

            outCfunction(
                outfile=os.path.join(Ccodesdir,name+".c"),
                includes=includes,
                desc=desc,
                c_type=c_type, name=name, params=params,
                enableCparameters=False,
                body=body)


        desc = "Add extrinsic curvature source term for tau_tilde"
        name = "calculate_tau_tilde_source_term_extrinsic_curv"    

        prestring += "const double KDD00 = curv->Kxx;\n"
        prestring += "const double KDD01 = curv->Kxy;\n"
        prestring += "const double KDD02 = curv->Kxz;\n"

        prestring += "const double KDD11 = curv->Kyy;\n"
        prestring += "const double KDD12 = curv->Kyz;\n"

        prestring += "const double KDD22 = curv->Kzz;\n"

        loopstring = prestring

        body = outputC(vars_rhs[-1], vars_to_write[-1], params=outCparams, 
                       filename="returnstring", prestring=loopstring)

        params  = "const primitive_quantities *restrict prims, "
        params  += "struct eos_parameters const *restrict eos, "
        params  += "const metric_quantities *restrict metric, "
        params  += "const extrinsic_curvature *restrict curv, "
        params  += "conservative_quantities *restrict cons"

        outCfunction(
            outfile=os.path.join(Ccodesdir,name+".c"),
            includes=includes,
            desc=desc,
            c_type=c_type, name=name, params=params,
            enableCparameters=False,
            body=body)


    
# GRHayL structs:

# typedef struct primitive_quantities {
# double rho, press, eps;
# double vx, vy, vz;
# double Bx, By, Bz;
# double entropy, Y_e, temperature;
# double enthalpy;
# } primitive_quantities;

# typedef struct conservative_quantities {
#   double rho, tau, Y_e;
#   double S_x, S_y, S_z;
#   double entropy;
# } conservative_quantities;

# typedef struct metric_quantities {
#   double adm_gxx, adm_gxy, adm_gxz;
#   double adm_gyy, adm_gyz, adm_gzz;
#   double adm_gupxx, adm_gupxy, adm_gupxz;
#   double adm_gupyy, adm_gupyz, adm_gupzz;
#   double betax, betay, betaz;
#   double lapse, lapseinv;
#   double psi2, psi4, psi6;
#   double psi4inv, lapseinv2;
#   double g4dn[4][4],g4up[4][4];
# } metric_quantities;

# typedef struct metric_derivatives {
#   double lapse[3];
#   double betax[3];
#   double betay[3];
#   double betaz[3];
#   double adm_gxx[3];
#   double adm_gxy[3];
#   double adm_gxz[3];
#   double adm_gyy[3];
#   double adm_gyz[3];
#   double adm_gzz[3];
# } metric_derivatives;