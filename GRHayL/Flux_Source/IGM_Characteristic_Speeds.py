# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import os
import sys

import GRMHD_equations_new_version as GRMHD  # NRPy+: Generate general relativistic magnetohydrodynamics equations

nrpy_dir_path = os.path.join("nrpy_old")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)

from sympy import symbols

from nrpy.c_codegen import c_codegen
from nrpy.c_function import CFunction
from nrpy.reference_metric import reference_metric
from nrpy.equations.general_relativity.ADM_to_BSSN import ADM_to_BSSN
from nrpy.equations.general_relativity.BSSN_quantities import BSSNQuantities
from nrpy.equations.grhd.characteristic_speeds import find_cmax_cmin
from nrpy.indexedexp import declarerank1, declarerank2, zerorank2
from nrpy.params import set_parval_from_str

# Set coordinate system to dst_basis
# set_parval_from_str("reference_metric::CoordSystem", "Cartesian")


def Cfunction__GRMHD_characteristic_speeds(
    Ccodesdir,
    includes=None,
    formalism="ADM",
):

    sqrt4pi = symbols("SQRT_4_PI")
    alpha_face = symbols("alpha_face")
    if formalism == "BSSN":
        cf_face = symbols("cf_face")
        h_faceDD = declarerank2("h_faceDD", symmetry="sym01", DIM=3)
        vet_faceU = declarerank1("vet_faceU", DIM=3)
        beta_faceU = vet_faceU

        AitoB = ADM_to_BSSN(None, None, None, None)
        Bq = BSSNQuantities()

        gamma_faceDD = zerorank2()
        for i in range(3):
            for j in range(3):
                gamma_faceDD[i][j] = (
                    AitoB.gammaDD[i][j]
                    .subs(Bq.hDD[i][j], h_faceDD[i][j])
                    .subs(Bq.cf, cf_face)
                )

    else:
        beta_faceU = declarerank1("beta_faceU")
        gamma_faceDD = declarerank2("gamma_faceDD", symmetry="sym01")

    # We'll need some more gridfunctions, now, to represent the reconstructions of BU and ValenciavU
    # on the right and left faces
    u4rU = declarerank1("u4rU", DIM=4)
    u4lU = declarerank1("u4lU", DIM=4)

    BrU = declarerank1("BrU")
    BlU = declarerank1("BlU")

    rho_b_r = symbols("rhob_r")
    rho_b_l = symbols("rhob_l")

    h_r = symbols("h_r")
    h_l = symbols("h_l")

    cs2_r = symbols("cs2_r")
    cs2_l = symbols("cs2_l")

    GRMHD.compute_smallb4U(gamma_faceDD, beta_faceU, alpha_face, u4rU, BrU, sqrt4pi)
    GRMHD.compute_smallbsquared(gamma_faceDD, beta_faceU, alpha_face, GRMHD.smallb4U)
    smallbsquared_r = GRMHD.smallbsquared

    GRMHD.compute_smallb4U(gamma_faceDD, beta_faceU, alpha_face, u4lU, BlU, sqrt4pi)
    GRMHD.compute_smallbsquared(gamma_faceDD, beta_faceU, alpha_face, GRMHD.smallb4U)
    smallbsquared_l = GRMHD.smallbsquared

    cmins_rhs = []
    cmaxs_rhs = []

    for flux_dirn in range(3):
        find_cmax_cmin(
            flux_dirn,
            gamma_faceDD,
            beta_faceU,
            alpha_face,
            smallbsquared_r,
            smallbsquared_l,
            u4rU,
            u4lU,
            rho_b_r,
            rho_b_l,
            h_r,
            h_l,
            cs2_r,
            cs2_l,
        )

        cmins_rhs.append(cmin)
        cmaxs_rhs.append(cmax)

    prims_NRPy_r = [
        "u4rU0",
        "u4rU1",
        "u4rU2",
        "u4rU3",
        "BrU0",
        "BrU1",
        "BrU2",
        "rhob_r",
    ]
    prims_NRPy_l = [
        "u4lU0",
        "u4lU1",
        "u4lU2",
        "u4lU3",
        "BlU0",
        "BlU1",
        "BlU2",
        "rhob_l",
    ]

    prims_GRHayL = ["u0", "vU[0]", "vU[1]", "vU[2]", "BU[0]", "BU[1]", "BU[2]", "rho"]

    prestring = r"""
double h_r, h_l, cs2_r, cs2_l;

ghl_compute_h_and_cs2(eos, prims_r, &h_r, &cs2_r);
ghl_compute_h_and_cs2(eos, prims_l, &h_l, &cs2_l);
"""

    for i in range(len(prims_NRPy_r)):
        if "v" in prims_GRHayL[i]:
            vel_r_str = prims_NRPy_r[0]
            vel_l_str = prims_NRPy_l[0]
            prestring += (
                "const double "
                + prims_NRPy_r[i]
                + " = prims_r->"
                + prims_GRHayL[i]
                + "*"
                + vel_r_str
                + ";\n"
            )
            prestring += (
                "const double "
                + prims_NRPy_l[i]
                + " = prims_l->"
                + prims_GRHayL[i]
                + "*"
                + vel_l_str
                + ";\n"
            )
        else:
            prestring += (
                "const double "
                + prims_NRPy_r[i]
                + " = prims_r->"
                + prims_GRHayL[i]
                + ";\n"
            )
            prestring += (
                "const double "
                + prims_NRPy_l[i]
                + " = prims_l->"
                + prims_GRHayL[i]
                + ";\n"
            )

    prestring += "const double " + str(alpha_face) + " = metric_face->lapse;\n"

    checker = []

    if formalism == "BSSN":
        prestring += (
            "const double "
            + str(cf_face)
            + " = metric_face_quantities->"
            + str(cf_face)
            + ";\n"
        )

        for i in range(3):
            vetU_var = vet_faceU[i]
            prestring += (
                "const double "
                + str(vetU_var)
                + " = metric_face_quantities->"
                + str(vetU_var)
                + ";\n"
            )

        for i in range(3):
            for j in range(3):
                hDD_var = h_faceDD[i][j]
                if hDD_var in checker:
                    continue
                prestring += (
                    "const double "
                    + str(hDD_var)
                    + " = metric_face_quantities->"
                    + str(hDD_var)
                    + ";\n"
                )
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

    cmins = ["cmin_dirn0", "cmin_dirn1", "cmin_dirn2"]

    cmaxs = ["cmax_dirn0", "cmax_dirn1", "cmax_dirn2"]

    c_type = "void"

    params = "const primitive_quantities *restrict prims_r, "
    params += "const primitive_quantities *restrict prims_l, "
    params += "const eos_parameters *restrict eos, "
    params += "const metric_quantities *restrict metric_face, "

    for flux_dirn in range(3):
        cmin_param = "double *" + cmins[flux_dirn] + ", "
        cmax_param = "double *" + cmaxs[flux_dirn]

        write_speeds_str = ["*" + cmins[flux_dirn], "*" + cmaxs[flux_dirn]]
        write_speeds_rhs_str = [cmins_rhs[flux_dirn], cmaxs_rhs[flux_dirn]]

        body = c_codegen(
            write_speeds_rhs_str,
            write_speeds_str,
            verbose=False,
            prestring=prestring,
            cse_sorting=True,
        )

        desc = "Compute the characteristic speeds in direction " + str(flux_dirn)
        name = "ghl_calculate_characteristic_speed_dirn" + str(flux_dirn)

        cfunc = CFunction(
            cfunc_type=c_type,
            includes=includes,
            desc=desc,
            name=name,
            params=params + cmin_param + cmax_param,
            body=body,
        )
        with open(os.path.join(Ccodesdir, name + ".c")) as f:
            f.write(cfunc.generate_full_function()[2])
