from typing import Tuple, List
import sympy as sp
import nrpy.indexedexp as ixp


# This calculates the stress-energy tensor
# for a given set of primitive variables
# using the equation
# T^{\mu\nu} = (rho h + b^2) u^\mu u^\nu
#            + (P + b^2/2) g^{\mu\nu}
#            - b^\mu b^nu
def calculate_T4UU(
    g4UU: List[List[sp.Expr]],
    rho: sp.Expr,
    press: sp.Expr,
    u4U: List[sp.Expr],
    smallb4U: List[sp.Expr],
    smallb_sq: sp.Expr,
    enthalpy: sp.Expr,
) -> List[List[sp.Expr]]:

    T4UU = ixp.zerorank2(dimension=4)
    for i in range(4):
        for j in range(4):
            T4UU[i][j] = (rho * enthalpy + smallb_sq) * u4U[i] * u4U[j] \
                       + (press + smallb_sq/2.0) * g4UU[i][j]           \
                       - smallb4U[i] * smallb4U[j]
    return T4UU

# alpha, alphainv, betaU, sqrt_detgamma, sqrt_detg,
# gammaDD, g4DD, g4UU
def initialize_ADM_metric_struct(
        metric_struct_name: str
    ) -> Tuple[sp.Expr, sp.Expr, List[sp.Expr], sp.Expr, sp.Expr, \
               List[List[sp.Expr]], List[List[sp.Expr]], List[List[sp.Expr]]]:

    # Step 1: set up space-time quantities
    alpha = sp.symbols(metric_struct_name + "->lapse")
    alphainv = sp.symbols(metric_struct_name + "->lapseinv")
    betaU = ixp.declarerank1("betaU")
    gammaDD = ixp.declarerank2("gammaDD", symmetry="sym01")
    sqrt_detgamma = sp.symbols(metric_struct_name + "->sqrt_detgamma")
    sqrt_detg = sp.symbols("sqrt_detg")
    g4DD = ixp.declarerank2("g4DD", symmetry="sym01", dimension=4)
    g4UU = ixp.declarerank2("g4UU", symmetry="sym01", dimension=4)

    return alpha, alphainv, betaU, sqrt_detgamma, sqrt_detg, \
           gammaDD, g4DD, g4UU

# rho, press, Ye, ent, u0, vU, BU, u4U
def initialize_primitives_struct(
        prims_struct_name: str
    ) -> Tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr, sp.Expr, sp.Expr, \
               List[sp.Expr], List[sp.Expr], List[sp.Expr]]:

    rho   = sp.symbols(prims_struct_name + "->rho")
    press = sp.symbols(prims_struct_name + "->press")
    Ye    = sp.symbols(prims_struct_name + "->Y_e")
    ent   = sp.symbols(prims_struct_name + "->ent")
    u0    = sp.symbols(prims_struct_name + "->u0")
    vU    = ixp.zerorank1()
    BU    = ixp.zerorank1()
    
    u4U = ixp.zerorank1(dimension=4)
    u4U[0] = u0
    for i in range(1,4):
        vU[i-1] = sp.symbols(prims_struct_name + "->vU[" + str(i-1) + "]")
        BU[i-1] = sp.symbols(prims_struct_name + "->BU[" + str(i-1) + "]")
        u4U[i] = vU[i-1] * u0

    return rho, press, Ye, ent, u0, vU, BU, u4U

# b4U, b4D, b^2
def initialize_smallb_vars(
        alphainv: sp.Expr,
        g4DD: List[List[sp.Expr]],
        u4U: List[sp.Expr],
        u4D: List[sp.Expr],
        BU: List[sp.Expr]
    ) -> Tuple[List[sp.Expr], List[sp.Expr], sp.Expr]:

    smallb4U = ixp.zerorank1(dimension=4)
    smallb4D = ixp.zerorank1(dimension=4)
    smallb_sq = sp.sympify(0)
    for i in range(3):
        smallb4U[0] += u4D[i+1] * BU[i] * alphainv

    for i in range(1,4):
        smallb4U[i] = (BU[i-1] * alphainv + smallb4U[0] * u4U[i]) / u4U[0]

    for i in range(4):
        for j in range(4):
            smallb4D[i] += g4DD[i][j] * smallb4U[j]

    for i in range(4):
        smallb_sq += smallb4D[i] * smallb4U[j]

    return smallb4U, smallb4D, smallb_sq

# rho_*, tau, SD, Ye_*, ent_*
def compute_conservatives(
        alpha: sp.Expr,
        sqrt_detg: sp.Expr,
        g4DD: List[List[sp.Expr]],
        rho: sp.Expr,
        u4U: List[sp.Expr],
        Ye: sp.Expr,
        ent: sp.Expr,
        T4UU: List[List[sp.Expr]]
    ) -> Tuple[sp.Expr, sp.Expr, List[sp.Expr], sp.Expr, sp.Expr]:

    rho_star = sqrt_detg*rho*u4U[0]
    tau = alpha*sqrt_detg*T4UU[0][0] - rho_star
    Ye_star = Ye*rho_star
    ent_star = sqrt_detg*ent*u4U[0]

    T4UD = ixp.zerorank2(dimension=4)
    for i in range(4):
        for j in range(4):
            for k in range(4):
                T4UD[i][j] += g4DD[j][k] * T4UU[i][k]

    SD = ixp.zerorank1()
    for i in range(3):
        SD[i] = sqrt_detg*T4UD[0][i+1]

    return rho_star, tau, SD, Ye_star, ent_star
