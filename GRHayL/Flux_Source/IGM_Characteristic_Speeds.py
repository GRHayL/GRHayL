# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import os,sys
import GRMHD_equations_new_version as GRMHD    # NRPy+: Generate general relativistic magnetohydrodynamics equations

nrpy_dir_path = os.path.join("nrpy/")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)

from outputC import outputC, outCfunction # NRPy+: Core C code output module
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import reference_metric as rfm   # NRPy+: Reference metric support

thismodule = __name__

CoordSystem = "Cartesian"

# Set coordinate system to dst_basis
par.set_parval_from_str("reference_metric::CoordSystem", CoordSystem)
rfm.reference_metric()

# We'll write this as a function so that we can calculate the expressions on-demand for any choice of i
def find_cp_cm(flux_dirn, g4UU,
               smallbsquared, u4U, rho_b, h, cs2):
    # Outputs: cplus,cminus
    vA2 = smallbsquared / (rho_b*h + smallbsquared)
    v02 = vA2 + (cs2)*(1 - vA2)

    a = (1 - v02)*(u4U[0]**2) - v02*g4UU[0][0]
    b = 2*v02*g4UU[flux_dirn+1][0] - 2*u4U[flux_dirn+1]*u4U[0]*(1 - v02)
    c = (1 - v02)*(u4U[flux_dirn+1]**2) - v02*g4UU[flux_dirn+1][flux_dirn+1]

    # Now, we are free to solve the quadratic equation as usual. We take care to avoid passing a
    # negative value to the sqrt function.
    detm = b*b - sp.sympify(4)*a*c

    import Min_Max_and_Piecewise_Expressions as noif
    detm = sp.sqrt(noif.max_noif(sp.sympify(0),detm))
    global cplus,cminus
    # note that these correspond to a single interface, left or right
    cplus_tmp  = sp.Rational(1,2)* (-b/a + detm/a)
    cminus_tmp = sp.Rational(1,2)*-( b/a + detm/a)

    cminus = noif.min_noif(cplus_tmp, cminus_tmp)
    cplus  = noif.max_noif(cplus_tmp, cminus_tmp)

    # the above in C code
    # if (cplus < cminus) {
    # CCTK_REAL cp = cminus;
    # cminus = cplus;
    # cplus = cp;


# We'll write this as a function, and call it within HLLE_solver, below.
def find_cmax_cmin(flux_dirn, gamma_faceDD, beta_faceU, alpha_face,
                   smallbsquared_r, smallbsquared_l,
                   u4U_r, u4U_l, rho_b_r, rho_b_l,
                   h_r, h_l, cs2_r, cs2_l):
    # Inputs:  flux direction flux_dirn, Inverse metric g4_faceUU, shift beta_faceU,
    #          lapse alpha_face, metric determinant gammadet_face
    # Outputs: maximum and minimum characteristic speeds cmax and cmin
    # First, we need to find the characteristic speeds on each face
    import BSSN.ADMBSSN_tofrom_4metric as AB4m
    AB4m.g4UU_ito_BSSN_or_ADM("ADM",gamma_faceDD, beta_faceU, alpha_face)

    # Original needed for GRMHD
    find_cp_cm(flux_dirn, AB4m.g4UU, smallbsquared_r, u4U_r, rho_b_r, h_r, cs2_r)
    cpr = cplus
    cmr = cminus

    find_cp_cm(flux_dirn, AB4m.g4UU, smallbsquared_l, u4U_l, rho_b_l, h_l, cs2_l)
    cpl = cplus
    cml = cminus

    # The following algorithms have been verified with random floats:

    global cmax,cmin
#   // Then compute cmax, cmin. This is required for the HLL flux.
#   original C code
#   CCTK_REAL cmaxL =  MAX(0.0,MAX(cplusl,cplusr));
#   CCTK_REAL cminL = -MIN(0.0,MIN(cminusl,cminusr));

    # Now, we need to set cmax to the larger of cpr,cpl, and 0
    import Min_Max_and_Piecewise_Expressions as noif
    cmax = noif.max_noif(0.0, noif.max_noif(cpl, cpr))
    # And then, set cmin to the smaller of cmr,cml, and 0
    cmin = -noif.min_noif(0.0, noif.min_noif(cml, cmr))

def Cfunction__GRMHD_characteristic_speeds(Ccodesdir, includes=None, formalism="ADM", outCparams = "outCverbose=False,CSE_sorting=True"):

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

    rho_b_r = sp.symbols("rhob_r")
    rho_b_l = sp.symbols("rhob_l")

    h_r = sp.symbols("h_r")
    h_l = sp.symbols("h_l")

    cs2_r = sp.symbols("cs2_r")
    cs2_l = sp.symbols("cs2_l")

    GRMHD.compute_smallb4U(gamma_faceDD, beta_faceU, alpha_face, u4rU, BrU, sqrt4pi)
    GRMHD.compute_smallbsquared(gamma_faceDD, beta_faceU, alpha_face, GRMHD.smallb4U)
    smallbsquared_r = GRMHD.smallbsquared

    GRMHD.compute_smallb4U(gamma_faceDD, beta_faceU, alpha_face, u4lU, BlU, sqrt4pi)
    GRMHD.compute_smallbsquared(gamma_faceDD, beta_faceU, alpha_face, GRMHD.smallb4U)
    smallbsquared_l = GRMHD.smallbsquared


    cmins_rhs = []
    cmaxs_rhs = []

    for flux_dirn in range(3):
        find_cmax_cmin(flux_dirn, gamma_faceDD, beta_faceU, alpha_face,
                       smallbsquared_r, smallbsquared_l,
                       u4rU, u4lU, rho_b_r, rho_b_l,
                       h_r, h_l, cs2_r, cs2_l,)

        cmins_rhs.append(cmin)
        cmaxs_rhs.append(cmax)


    prims_NRPy_r = ["u4rU0", "u4rU1", "u4rU2", "u4rU3", "BrU0", "BrU1", "BrU2", "rhob_r"]
    prims_NRPy_l = ["u4lU0", "u4lU1", "u4lU2", "u4lU3", "BlU0", "BlU1", "BlU2", "rhob_l"]

    prims_GRHayL = ["u0", "vU[0]", "vU[1]", "vU[2]", "BU[0]", "BU[1]", "BU[2]", "rho"]

    prestring = r"""
double h_r, h_l, cs2_r, cs2_l;

eos->compute_h_and_cs2(eos, prims_r, &h_r, &cs2_r);
eos->compute_h_and_cs2(eos, prims_l, &h_l, &cs2_l);
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


    cmins = ["cmin_dirn0",
             "cmin_dirn1",
             "cmin_dirn2"]

    cmaxs = ["cmax_dirn0",
             "cmax_dirn1",
             "cmax_dirn2"]

    c_type = "void"

    params  =  "const primitive_quantities *restrict prims_r, "
    params  += "const primitive_quantities *restrict prims_l, "
    params  += "struct eos_parameters const *restrict eos, "
    params  += "const metric_quantities *restrict metric_face, "


    for flux_dirn in range(3):
        cmin_param = "double *"+cmins[flux_dirn]+", "
        cmax_param = "double *"+cmaxs[flux_dirn]

        write_speeds_str = ["*"+cmins[flux_dirn], "*"+cmaxs[flux_dirn]]
        write_speeds_rhs_str = [cmins_rhs[flux_dirn], cmaxs_rhs[flux_dirn]]

        body = outputC(write_speeds_rhs_str, write_speeds_str, params=outCparams,
                       filename="returnstring", prestring=prestring)

        desc = "Compute the characteristic speeds in direction " + str(flux_dirn)
        name = "ghl_calculate_characteristic_speed_dirn" + str(flux_dirn)

        outCfunction(
            outfile=os.path.join(Ccodesdir,name+".c"),
            includes=includes,
            desc=desc,
            name=name,
            params=params + cmin_param + cmax_param,
            body= body,
            enableCparameters=False)
