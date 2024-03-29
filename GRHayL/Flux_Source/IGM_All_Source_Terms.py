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
    prims_GRHayL = ["u0", "vU[0]*u4U0", "vU[1]*u4U0", "vU[2]*u4U0", "BU[0]", "BU[1]", "BU[2]", "press", "rho"]

    prestring = r"""
double h, cs2;

ghl_compute_h_and_cs2(eos, prims, &h, &cs2);
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


        prestring += "const double betaU0 = metric->betaU[0];\n"
        prestring += "const double betaU1 = metric->betaU[1];\n"
        prestring += "const double betaU2 = metric->betaU[2];\n"

        prestring += "const double gammaDD00 = metric->gammaDD[0][0];\n"
        prestring += "const double gammaDD01 = metric->gammaDD[0][1];\n"
        prestring += "const double gammaDD02 = metric->gammaDD[0][2];\n"

        prestring += "const double gammaDD11 = metric->gammaDD[1][1];\n"
        prestring += "const double gammaDD12 = metric->gammaDD[1][2];\n"

        prestring += "const double gammaDD22 = metric->gammaDD[2][2];\n"

        prestring += "const double KDD00 = curv->K[0][0];\n"
        prestring += "const double KDD01 = curv->K[0][1];\n"
        prestring += "const double KDD02 = curv->K[0][2];\n"

        prestring += "const double KDD11 = curv->K[1][1];\n"
        prestring += "const double KDD12 = curv->K[1][2];\n"

        prestring += "const double KDD22 = curv->K[2][2];\n"

        for i in range(3):
            prestring += f"const double alpha_dD{i}     = metric_derivs_{chr(ord('x')+i)}->lapse;\n"
            prestring += f"const double betaU_dD0{i}    = metric_derivs_{chr(ord('x')+i)}->betaU[0];\n"
            prestring += f"const double betaU_dD1{i}    = metric_derivs_{chr(ord('x')+i)}->betaU[1];\n"
            prestring += f"const double betaU_dD2{i}    = metric_derivs_{chr(ord('x')+i)}->betaU[2];\n"
            prestring += f"const double gammaDD_dD00{i} = metric_derivs_{chr(ord('x')+i)}->gammaDD[0][0];\n"
            prestring += f"const double gammaDD_dD01{i} = metric_derivs_{chr(ord('x')+i)}->gammaDD[0][1];\n"
            prestring += f"const double gammaDD_dD02{i} = metric_derivs_{chr(ord('x')+i)}->gammaDD[0][2];\n"
            prestring += f"const double gammaDD_dD11{i} = metric_derivs_{chr(ord('x')+i)}->gammaDD[1][1];\n"
            prestring += f"const double gammaDD_dD12{i} = metric_derivs_{chr(ord('x')+i)}->gammaDD[1][2];\n"
            prestring += f"const double gammaDD_dD22{i} = metric_derivs_{chr(ord('x')+i)}->gammaDD[2][2];\n"

        desc     = f"Add source terms for Stilde and tau_tilde"
        c_type   = "void"
        name     = f"ghl_calculate_source_terms"
        params   = "const ghl_eos_parameters *restrict eos, "
        params  += "ghl_primitive_quantities *restrict prims, "
        params  += "const ghl_metric_quantities *restrict metric, "
        params  += "const ghl_metric_quantities *restrict metric_derivs_x, "
        params  += "const ghl_metric_quantities *restrict metric_derivs_y, "
        params  += "const ghl_metric_quantities *restrict metric_derivs_z, "
        params  += "const ghl_extrinsic_curvature *restrict curv, "
        params  += "ghl_conservative_quantities *restrict cons"

        vars_to_write = ["cons->SD[0]", "cons->SD[1]", "cons->SD[2]", "cons->tau"]

        tau_tilde_source = GRMHD.tau_tilde_source_term_2_3_A + \
                           GRMHD.tau_tilde_source_term_2_3_B + \
                           GRMHD.tau_tilde_source_term_2_3_C + \
                           GRMHD.tau_tilde_source_term_1_3

        vars_rhs = [GRMHD.S_tilde_source_termD[0],
                    GRMHD.S_tilde_source_termD[1],
                    GRMHD.S_tilde_source_termD[2],
                    tau_tilde_source]

        body = outputC(vars_rhs, vars_to_write,
                       params=outCparams,
                       filename="returnstring", prestring=prestring)

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
#   double g4DD[4][4],g4UU[4][4];
# } metric_quantities;
