from typing import Tuple, List
import sympy as sp
import nrpy.indexedexp as ixp
import ghl_helpers as ghl

# This calculates the HLL flux for a variable
# U using the formula
# F = (c_min F_r + c_max F_l - c_min c_max (U_r - U_l) ) / (c_min + c_max)
def calculate_HLL_flux(
    cmin: sp.Expr,
    cmax: sp.Expr,
    flux_r: sp.Expr,
    var_r: sp.Expr,
    flux_l: sp.Expr,
    var_l: sp.Expr
) -> sp.Expr:

    return (cmin*flux_r + cmax*flux_l - cmin*cmax*(var_r - var_l) ) / (cmin + cmax)

def calculate_fluxes_rhs() -> Tuple[sp.Expr, sp.Expr, List[sp.Expr]]:

    # Step 1: set up space-time quantities
    alpha, alphainv, betaU, sqrt_detgamma, sqrt_detg, gammaDD, g4DD, g4UU = \
        ghl.initialize_ADM_metric_struct("ADM_metric")

    # Step 2: set up fluid quantities
    cmin = sp.symbols("cmin")
    cmax = sp.symbols("cmax")
    
    # Step 2.1: initialize variables for the right side
    rho_r, press_r, Ye_r, ent_r, u0_r, vU_r, BU_r, u4U_r = \
        ghl.initialize_primitives_struct("prims_r")

    h_r     = sp.symbols("h_r")
    u4D_r = ixp.zerorank1(dimension=4)
    for i in range(4):
        for j in range(4):
            u4D_r[i] += g4DD[i][j] * u4U_r[j]

    smallb4U_r, smallb4D_r, smallb_sq_r = \
        ghl.initialize_smallb_vars(alphainv, g4DD, u4U_r, u4D_r, BU_r)

    # Step 2.2: initialize variables for the left side
    rho_l, press_l, Ye_l, ent_l, u0_l, vU_l, BU_l, u4U_l = \
    ghl.initialize_primitives_struct("prims_l")

    h_l     = sp.symbols("h_l")
    u4D_l = ixp.zerorank1(dimension=4)
    for i in range(4):
        for j in range(4):
            u4D_l[i] += g4DD[i][j] * u4U_l[j]

    smallb4U_l, smallb4D_l, smallb_sq_l = \
        ghl.initialize_smallb_vars(alphainv, g4DD, u4U_l, u4D_l, BU_l)

    # Step 3: Compute Tmunu
    T4UU_r = ghl.calculate_T4UU(g4UU, rho_r, press_r, u4U_r, smallb4U_r, smallb_sq_r, h_r)
    T4UU_l = ghl.calculate_T4UU(g4UU, rho_l, press_l, u4U_l, smallb4U_l, smallb_sq_l, h_l)

    # Step 3: Define directional components
    #   We do this explicitly so that we don't have to generate
    #   different functions for each flux direction with the
    #   understanding that there is a variable "direction" in
    #   the C code that will access the correct array element.

    # directional 4-metric
    gdir_symbol_0 = sp.symbols("metric_aux.g4UU[0][direction+1]")
    gdir_symbol_1 = sp.symbols("metric_aux.g4UU[1][direction+1]")
    gdir_symbol_2 = sp.symbols("metric_aux.g4UU[2][direction+1]")
    gdir_symbol_3 = sp.symbols("metric_aux.g4UU[3][direction+1]")
    g4UU_dir = ixp.zerorank1(dimension=4)
    g4UU_dir[0] = gdir_symbol_0
    g4UU_dir[1] = gdir_symbol_1
    g4UU_dir[2] = gdir_symbol_2
    g4UU_dir[3] = gdir_symbol_3

    # This variable = g^{d \mu} g_{\mu \nu} where d is the
    # direction of the flux
    g4UD_dir = ixp.zerorank1(dimension=4)
    for i in range(4):
        for j in range(4):
            g4UD_dir[i] += g4UU_dir[j] * g4DD[i][j]

    # directional 3-velocity
    v_dir_r = sp.symbols("prims_r->vU[direction]")
    v_dir_l = sp.symbols("prims_l->vU[direction]")

    # directional B and smallb
    B_dir_r = sp.symbols("prims_r->BU[direction]")
    B_dir_l = sp.symbols("prims_l->BU[direction]")
    smallb_dir_r = B_dir_r * alphainv / u4U_r[0] + smallb4U_r[0] * v_dir_r
    smallb_dir_l = B_dir_l * alphainv / u4U_l[0] + smallb4U_l[0] * v_dir_l

    # Directional stress-energy vector T^{0i}
    # Note that we have converted u^i -> u^0 v^i to use the directional variable
    T_0i_dir_r = (rho_r * h_r + smallb_sq_r) * u4U_r[0] * u4U_r[0] * v_dir_r \
               + (press_r + smallb_sq_r/2.0) * g4UU_dir[0]                          \
               - smallb4U_r[0] * smallb_dir_r
    T_0i_dir_l = (rho_l * h_l + smallb_sq_l) * u4U_l[0] * u4U_l[0] * v_dir_l \
               + (press_l + smallb_sq_l/2.0) * g4UU_dir[0]                          \
               - smallb4U_l[0] * smallb_dir_l

    # directional stress-energy tensor T^i_j
    T_ij_dir_r = ixp.zerorank1()
    T_ij_dir_l = ixp.zerorank1()
    for i in range(3):
        T_ij_dir_r[i] = (rho_r * h_r + smallb_sq_r) * u4U_r[0] * v_dir_r  * u4D_l[i] \
                      + (press_r + smallb_sq_r/2.0) * g4UD_dir[i+1]                         \
                      - smallb_dir_r * smallb4D_r[i]
        T_ij_dir_l[i] = (rho_l * h_l + smallb_sq_l) * u4U_l[0] * v_dir_l  * u4D_l[i] \
                      + (press_l + smallb_sq_l/2.0) * g4UD_dir[i+1]                         \
                      - smallb_dir_l * smallb4D_l[i]

    # Step 4: compute conservative variables
    rho_star_r, tau_r, SD_r, Ye_star_r, ent_star_r = \
        ghl.compute_conservatives(alpha, sqrt_detg, g4DD, rho_r, u4U_r, Ye_r, ent_r, T4UU_r)

    rho_star_l, tau_l, SD_l, Ye_star_l, ent_star_l = \
        ghl.compute_conservatives(alpha, sqrt_detg, g4DD, rho_l, u4U_l, Ye_l, ent_l, T4UU_l)

    # Step 5: compute fluxes on each side
    rho_flux_r = rho_star_r * v_dir_r
    rho_flux_l = rho_star_l * v_dir_l

    tau_flux_r = alpha*sqrt_detg*T_0i_dir_r - rho_flux_r
    tau_flux_l = alpha*sqrt_detg*T_0i_dir_l - rho_flux_l

    S_flux_r = ixp.zerorank1()
    S_flux_l = ixp.zerorank1()
    for i in range(3):
        S_flux_r[i] = sqrt_detg*T_ij_dir_r[i]
        S_flux_l[i] = sqrt_detg*T_ij_dir_l[i]

    Ye_flux_r = rho_flux_r * Ye_r
    Ye_flux_l = rho_flux_l * Ye_l

    ent_flux_r = ent_star_r * v_dir_r
    ent_flux_l = ent_star_l * v_dir_l

    # Step 6: compute rhs contributions due to fluxes
    rho_flux_rhs = calculate_HLL_flux( \
        cmin, cmax, rho_flux_r, rho_star_r, rho_flux_l, rho_star_l)

    tau_flux_rhs = calculate_HLL_flux( \
        cmin, cmax, tau_flux_r, tau_r, tau_flux_l, tau_l)

    S_flux_rhs = ixp.zerorank1()
    for i in range(3):
        S_flux_rhs[i] = calculate_HLL_flux( \
                cmin, cmax, S_flux_r[i], SD_r[i], S_flux_l[i], SD_l[i])

    Ye_flux_rhs = calculate_HLL_flux( \
        cmin, cmax, Ye_flux_r, Ye_star_r, Ye_flux_l, Ye_star_l)

    ent_flux_rhs = calculate_HLL_flux( \
        cmin, cmax, ent_flux_r, ent_star_r, ent_flux_l, ent_star_l)

    return rho_flux_rhs, tau_flux_rhs, S_flux_rhs, Ye_flux_rhs, ent_flux_rhs
