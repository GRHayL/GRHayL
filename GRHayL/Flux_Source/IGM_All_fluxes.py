# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import os,sys
import GRMHD_equations_new_version as GRMHD    # NRPy+: Generate general relativistic magnetohydrodynamics equations

nrpy_dir_path = os.path.join("nrpy/")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)

from outputC import outputC, outCfunction # NRPy+: Core C code output module
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import finite_difference as fin  # NRPy+: Finite difference C code generation module
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm   # NRPy+: Reference metric support

#Step 0: Set the spatial dimension parameter to 3.
par.set_parval_from_str("grid::DIM", 3)
DIM = par.parval_from_str("grid::DIM")

CoordSystem = "Cartesian"

# Set coordinate system to dst_basis
par.set_parval_from_str("reference_metric::CoordSystem", CoordSystem)
rfm.reference_metric()


# We'll rewrite this assuming that we've passed the entire reconstructed
# gridfunctions. You could also do this with only one point, but then you'd
# need to declare everything as a Cparam in NRPy+

def calculate_GRMHD_Tmunu_and_contractions(formalism, flux_dirn, mom_comp,
                                           gammaDD,betaU,alpha,
                                           rho_b, P, h, u4U, BU, S=None, Y_e=None):
    GRMHD.set_up_base_vars(formalism=formalism)

    GRMHD.compute_vU_from_u4U__no_speed_limit(u4U)

    # First compute stress-energy tensor T4UU and T4UD:
    GRMHD.compute_sqrtgammaDET(gammaDD)
    GRMHD.compute_smallb4U(gammaDD, betaU, alpha, u4U, BU, GRMHD.sqrt4pi)
    GRMHD.compute_smallbsquared(gammaDD, betaU, alpha, GRMHD.smallb4U)

    # First compute stress-energy tensor T4UU and T4UD:
    GRMHD.compute_T4UU(gammaDD,betaU,alpha, rho_b, P, h, u4U, GRMHD.smallb4U, GRMHD.smallbsquared)
    GRMHD.compute_T4UD(gammaDD,betaU,alpha, GRMHD.T4UU)

    # Compute conservative variables in terms of primitive variables
    GRMHD.compute_rho_star(alpha, GRMHD.sqrtgammaDET, rho_b,u4U)
    GRMHD.compute_tau_tilde(alpha, GRMHD.sqrtgammaDET, GRMHD.T4UU,GRMHD.rho_star)
    GRMHD.compute_S_tildeD(alpha, GRMHD.sqrtgammaDET, GRMHD.T4UD)

    # Next compute fluxes of conservative variables
    GRMHD.compute_rho_star_fluxU(GRMHD.VU,GRMHD.rho_star)
    GRMHD.compute_tau_tilde_fluxU(alpha, GRMHD.sqrtgammaDET, GRMHD.VU, GRMHD.T4UU, GRMHD.rho_star)
    GRMHD.compute_S_tilde_fluxUD (alpha, GRMHD.sqrtgammaDET, GRMHD.T4UD)

    global U_rho_star, F_rho_star
    U_rho_star = GRMHD.rho_star
    F_rho_star = GRMHD.rho_star_fluxU[flux_dirn]

    global U_tau_tilde, F_tau_tilde
    U_tau_tilde = GRMHD.tau_tilde
    F_tau_tilde = GRMHD.tau_tilde_fluxU[flux_dirn]

    global U_S_tilde, F_S_tilde
    U_S_tilde = GRMHD.S_tildeD[mom_comp]
    F_S_tilde = GRMHD.S_tilde_fluxUD[flux_dirn][mom_comp]

    global U_S_star, F_S_star
    U_S_star = None
    F_S_star = None
    if S is not None:
        GRMHD.compute_S_star(alpha, GRMHD.sqrtgammaDET, S, u4U)
        GRMHD.compute_S_star_fluxU(GRMHD.VU, GRMHD.S_star)
        U_S_star = GRMHD.S_star
        F_S_star = GRMHD.S_star_fluxU[flux_dirn]

    global U_Y_e_star, F_Y_e_star
    U_Y_e_star = None
    F_Y_e_star = None
    if Y_e is not None:
        GRMHD.compute_Y_e_star(alpha, GRMHD.sqrtgammaDET, Y_e, rho_b, u4U)
        GRMHD.compute_Y_e_star_fluxU(GRMHD.VU, GRMHD.Y_e_star)
        U_Y_e_star = GRMHD.Y_e_star
        F_Y_e_star = GRMHD.Y_e_star_fluxU[flux_dirn]

def HLLE_solver(cmax, cmin, Fr, Fl, Ur, Ul):
    # This solves the Riemann problem for the mom_comp component of the momentum
    # flux StildeD in the flux_dirn direction.

    # st_j_flux = (c_\min f_R + c_\max f_L - c_\min c_\max ( st_j_r - st_j_l )) / (c_\min + c_\max)
    return (cmin*Fr + cmax*Fl - cmin*cmax*(Ur-Ul) )/(cmax + cmin)




def calculate_HLLE_fluxes(formalism, flux_dirn, alpha_face, gamma_faceDD, beta_faceU,
                          u4rU, u4lU, BrU, BlU,
                          rho_b_r, rho_b_l,
                          P_r, P_l,
                          h_r, h_l,
                          cmin, cmax,
                          S_r=None, S_l=None,
                          Y_e_r=None, Y_e_l=None):

    global Stilde_flux_HLLED, rho_star_HLLE_flux, tau_tilde_HLLE_flux, S_star_HLLE_flux, Y_e_star_HLLE_flux
    Stilde_flux_HLLED = ixp.zerorank1()
#     rescaled_Stilde_fluxD = ixp.zerorank1()
    for mom_comp in range(3):
        calculate_GRMHD_Tmunu_and_contractions(formalism, flux_dirn, mom_comp,
                                               gamma_faceDD,beta_faceU,
                                               alpha_face,
                                               rho_b_r, P_r, h_r, u4rU,
                                               BrU, S=S_r, Y_e=Y_e_r)

        F_S_tilde_r = F_S_tilde
        U_S_tilde_r = U_S_tilde

        if mom_comp==0:
            U_rho_star_r = U_rho_star
            F_rho_star_r = F_rho_star

            U_tau_tilde_r = U_tau_tilde
            F_tau_tilde_r = F_tau_tilde

            if S_r is not None:
                U_S_star_r = U_S_star
                F_S_star_r = F_S_star

            if Y_e_r is not None:
                U_Y_e_star_r = U_Y_e_star
                F_Y_e_star_r = F_Y_e_star

        calculate_GRMHD_Tmunu_and_contractions(formalism, flux_dirn, mom_comp,
                                               gamma_faceDD,beta_faceU,
                                               alpha_face,
                                               rho_b_l, P_l, h_l, u4lU,
                                               BlU, S=S_l, Y_e=Y_e_l)

        F_S_tilde_l = F_S_tilde
        U_S_tilde_l = U_S_tilde

        if mom_comp==0:
            U_rho_star_l = U_rho_star
            F_rho_star_l = F_rho_star

            U_tau_tilde_l = U_tau_tilde
            F_tau_tilde_l = F_tau_tilde

            # now calculate HLLE derived fluxes, and rescale all of them
            rho_star_HLLE_flux = HLLE_solver(cmax, cmin,
                                      F_rho_star_r, F_rho_star_l,
                                      U_rho_star_r, U_rho_star_l)

            tau_tilde_HLLE_flux = HLLE_solver(cmax, cmin,
                                      F_tau_tilde_r, F_tau_tilde_l,
                                      U_tau_tilde_r, U_tau_tilde_l)

            if S_l is not None:
                U_S_star_l = U_S_star
                F_S_star_l = F_S_star
                S_star_HLLE_flux = HLLE_solver(cmax, cmin,
                                    F_S_star_r, F_S_star_l,
                                    U_S_star_r, U_S_star_l)

            if Y_e_l is not None:
                U_Y_e_star_l = U_Y_e_star
                F_Y_e_star_l = F_Y_e_star
                Y_e_star_HLLE_flux = HLLE_solver(cmax, cmin,
                                      F_Y_e_star_r, F_Y_e_star_l,
                                      U_Y_e_star_r, U_Y_e_star_l)

        # Rescale the flux term, to be FD
        Stilde_flux_HLLED[mom_comp] = HLLE_solver(cmax, cmin,
                                      F_S_tilde_r, F_S_tilde_l,
                                      U_S_tilde_r, U_S_tilde_l)


def Cfunction__GRMHD_fluxes(Ccodesdir, formalism="ADM", includes=None, tabulated=False, entropy=False,
                            outCparams = "outCverbose=False,CSE_sorting=True"):

    sqrt4pi = sp.symbols("SQRT_4_PI")
    alpha_face = sp.symbols("alpha_face")
    if formalism=="BSSN":
        cf_face = sp.symbols("cf_face")
        h_faceDD = ixp.declarerank2("h_faceDD","sym01",DIM=3)
        vet_faceU = ixp.declarerank1("vet_faceU", DIM=3)
        beta_faceU = vet_faceU

        import BSSN.ADM_in_terms_of_BSSN as AitoB
        AitoB.ADM_in_terms_of_BSSN()

        import BSSN.BSSN_quantities as Bq

        gamma_faceDD = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                gamma_faceDD[i][j] = AitoB.gammaDD[i][j].subs(Bq.hDD[i][j], h_faceDD[i][j]).subs(Bq.cf, cf_face)

    else:
        beta_faceU = ixp.declarerank1("beta_faceU")
        gamma_faceDD = ixp.declarerank2("gamma_faceDD","sym01")


    # We'll need some more gridfunctions, now, to represent the reconstructions of BU and ValenciavU
    # on the right and left faces
    u4rU = ixp.declarerank1("u4rU", DIM=4)
    u4lU = ixp.declarerank1("u4lU", DIM=4)

    BrU = ixp.declarerank1("BrU")
    BlU = ixp.declarerank1("BlU")

    cmin_dirn0 = sp.symbols("cmin_dirn0")
    cmin_dirn1 = sp.symbols("cmin_dirn1")
    cmin_dirn2 = sp.symbols("cmin_dirn2")

    cmax_dirn0 = sp.symbols("cmax_dirn0")
    cmax_dirn1 = sp.symbols("cmax_dirn1")
    cmax_dirn2 = sp.symbols("cmax_dirn2")

    cmins = [cmin_dirn0, cmin_dirn1, cmin_dirn2]
    cmaxs = [cmax_dirn0, cmax_dirn1, cmax_dirn2]

    h_r = sp.symbols("h_r")
    h_l = sp.symbols("h_l")

    P_r = sp.symbols("P_r")
    P_l = sp.symbols("P_l")

    rho_b_r = sp.symbols("rhob_r")
    rho_b_l = sp.symbols("rhob_l")

    prims_NRPy_r = ["u4rU0", "u4rU1", "u4rU2", "u4rU3", "BrU0", "BrU1", "BrU2", "P_r", "rhob_r"]
    prims_NRPy_l = ["u4lU0", "u4lU1", "u4lU2", "u4lU3", "BlU0", "BlU1", "BlU2", "P_l", "rhob_l"]

    prims_GRHayL = ["u0", "vU[0]", "vU[1]", "vU[2]", "BU[0]", "BU[1]", "BU[2]", "press", "rho"]

    S_r = None
    S_l = None
    if entropy:
        S_r = sp.symbols("S_r")
        S_l = sp.symbols("S_l")
        prims_NRPy_r += ["S_r"]
        prims_NRPy_l += ["S_l"]
        prims_GRHayL += ["entropy"]

    Y_e_r = None
    Y_e_l = None
    if tabulated:
        Y_e_r = sp.symbols("Y_e_r")
        Y_e_l = sp.symbols("Y_e_l")
        prims_NRPy_r += ["Y_e_r"]
        prims_NRPy_l += ["Y_e_l"]
        prims_GRHayL += ["Y_e"]

    prestring = r"""
double h_r, h_l, cs2_r, cs2_l;

ghl_compute_h_and_cs2(eos, prims_r, &h_r, &cs2_r);
ghl_compute_h_and_cs2(eos, prims_l, &h_l, &cs2_l);
"""

    for i in range(len(prims_NRPy_r)):
        if "v" in prims_GRHayL[i]:
            vel_r_str = prims_NRPy_r[0]
            vel_l_str = prims_NRPy_l[0]
            prestring += "const double "+prims_NRPy_r[i]+" = prims_r->"+prims_GRHayL[i]+"*"+vel_r_str+";\n"
            prestring += "const double "+prims_NRPy_l[i]+" = prims_l->"+prims_GRHayL[i]+"*"+vel_l_str+";\n"
        else:
            prestring += "const double "+prims_NRPy_r[i]+" = prims_r->"+prims_GRHayL[i]+";\n"
            prestring += "const double "+prims_NRPy_l[i]+" = prims_l->"+prims_GRHayL[i]+";\n"

    prestring += "const double "+str(alpha_face)+" = metric_face->lapse;\n"

    checker = []

    if formalism=="BSSN":
        prestring += "const double "+str(cf_face)+" = metric_face_quantities->"+str(cf_face)+";\n"

        for i in range(3):
            vetU_var = vet_faceU[i]
            prestring += "const double "+str(vetU_var)+" = metric_face_quantities->"+str(vetU_var)+";\n"

        for i in range(3):
            for j in range(3):
                hDD_var = h_faceDD[i][j]
                if hDD_var in checker:
                    continue
                prestring += "const double "+str(hDD_var)+" = metric_face_quantities->"+str(hDD_var)+";\n"
                checker.append(hDD_var)

    else:
        prestring += "const double beta_faceU0 = metric_face->betaU[0];\n"
        prestring += "const double beta_faceU1 = metric_face->betaU[1];\n"
        prestring += "const double beta_faceU2 = metric_face->betaU[2];\n"

        prestring += "const double gamma_faceDD00 = metric_face->gammaDD[0][0];\n"
        prestring += "const double gamma_faceDD01 = metric_face->gammaDD[0][1];\n"
        prestring += "const double gamma_faceDD02 = metric_face->gammaDD[0][2];\n"

        prestring += "const double gamma_faceDD11 = metric_face->gammaDD[1][1];\n"
        prestring += "const double gamma_faceDD12 = metric_face->gammaDD[1][2];\n"

        prestring += "const double gamma_faceDD22 = metric_face->gammaDD[2][2];\n"

    #     for i in range(3):
    #             betaU_var = beta_faceU[i]
    #             prestring += "const double "+str(betaU_var)+" = metric_face_quantities->"+str(betaU_var)+";\n"

    #     for i in range(3):
    #         for j in range(3):
    #             gammaDD_var = gamma_faceDD[i][j]
    #             if gammaDD_var in checker:
    #                 continue
    #             prestring += "const double "+str(gammaDD_var)+" = metric_face_quantities->"+str(gammaDD_var)+";\n"
    #             checker.append(gammaDD_var)


    vars_to_write = ["cons->SD[0]", "cons->SD[1]", "cons->SD[2]", "cons->rho", "cons->tau"]
    if entropy:
        vars_to_write += ["cons->entropy"]
    if tabulated:
        vars_to_write += ["cons->Y_e"]

    c_type = "void"

    params  =  "const primitive_quantities *restrict prims_r, "
    params  += "const primitive_quantities *restrict prims_l, "
    params  += "struct eos_parameters const *restrict eos, "
    params  += "const metric_quantities *restrict metric_face, "

    for flux_dirn in range(3):
        cmin_cmax_str = "const double "+str(cmins[flux_dirn])+", const double "+str(cmaxs[flux_dirn])+", "

    #     calc_char_speeds_params_str = "(prims_r, prims_l, eos, metric_face, &"+str(cmins[flux_dirn])+", &"+str(cmaxs[flux_dirn])+")"


    #     calc_char_speeds_func_str = "calculate_characteristic_speed_dirn" + str(flux_dirn)

        calculate_HLLE_fluxes(formalism, flux_dirn, alpha_face, gamma_faceDD, beta_faceU,
                              u4rU, u4lU, BrU, BlU,
                              rho_b_r, rho_b_l,
                              P_r, P_l,
                              h_r, h_l,
                              cmins[flux_dirn], cmaxs[flux_dirn],
                              S_r=S_r, S_l=S_l,
                              Y_e_r=Y_e_r, Y_e_l=Y_e_l)

        vars_rhs = [Stilde_flux_HLLED[0],
                    Stilde_flux_HLLED[1],
                    Stilde_flux_HLLED[2],
                    rho_star_HLLE_flux,
                    tau_tilde_HLLE_flux]
        if entropy:
            vars_rhs += [S_star_HLLE_flux]
        if tabulated:
            vars_rhs += [Y_e_star_HLLE_flux]

        body = outputC(vars_rhs, vars_to_write, params=outCparams,
                   filename="returnstring", prestring=prestring)

    #     prestring=(cmin_cmax_str+ calc_char_speeds_func_str+
    #                                calc_char_speeds_params_str+";\n\n"+
    #                                prestring)

        desc = "Compute the HLLE-derived fluxes on the left face in the " + str(flux_dirn) + "direction for all components."
        name = "ghl_calculate_HLLE_fluxes_dirn" + str(flux_dirn) + "_" + Ccodesdir

        outCfunction(
            outfile=os.path.join(Ccodesdir,name+".c"),
            includes=includes,
            desc=desc,
            name=name,
            params=params+cmin_cmax_str+"conservative_quantities *restrict cons",
            body= body,
            enableCparameters=False)
