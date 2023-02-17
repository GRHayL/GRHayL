# As documented in the NRPy+ tutorial module
#   Tutorial-GRFFE_Equations-Cartesian-new.ipynb
#   this module will construct useful quantities
#   for the GiRaFFE implementation of
#   general relativistic force-free
#   electrodynamics (GRFFE)
#
# Authors: Zachariah B. Etienne
#          zachetie **at** gmail **dot* com
#          Patrick Nelson
#          Terrence Pierre Jacques

# Step 1: Import needed core NRPy+ modules
# Import needed Python modules
# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import os,sys
nrpy_dir_path = os.path.join("nrpy/")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)

import NRPy_param_funcs as par   # NRPy+: Parameter interface
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm        # NRPy+: Reference metric support
import sympy as sp
from outputC import nrpyAbs      # NRPy+: Core C code output module

#Step 0: Set the spatial dimension parameter to 3.
par.set_parval_from_str("grid::DIM", 3)
DIM = par.parval_from_str("grid::DIM")

CoordSystem = "Cartesian"

# Set coordinate system to dst_basis
par.set_parval_from_str("reference_metric::CoordSystem", CoordSystem)
rfm.reference_metric()

thismodule = __name__
TINYDOUBLE = par.Cparameters("#define",thismodule,"TINYDOUBLE",1e-100)
# The default value isn't terribly important here, since we can overwrite in the main C code
GAMMA_SPEED_LIMIT = par.Cparameters("REAL", thismodule, "GAMMA_SPEED_LIMIT", 10.0)  # Default value based on IllinoisGRMHD.
# GiRaFFE default = 2000.0

# M_PI  = par.Cparameters("#define",thismodule,["M_PI"], "")
# sqrt4pi = par.Cparameters("REAL",thismodule,"sqrt4pi","sqrt(4.0*M_PI)")
sqrt4pi = sp.symbols("SQRT_4_PI")

import BSSN.BSSN_quantities as Bq
# par.set_parval_from_str("BSSN.BSSN_quantities::EvolvedConformalFactor_cf", "W")

def set_up_base_vars(formalism="BSSN"):
    global u4U, betaU, betaU_dD, alpha, alpha_dD, sqrt4pi, KDD, gammaDD, sqrtgammaDET, gammaDD_dD
    global BU, ValenciavU, VU, rho_b, P, h, epsilon 
    
    if formalism == "BSSN":
        Bq.BSSN_basic_tensors()
        Bq.gammabar__inverse_and_derivs()
        Bq.detgammabar_and_derivs()
        Bq.betaU_derivs()
        
        betaU    = Bq.betaU
        betaU_dD = Bq.betaU_dD

        import BSSN.ADM_in_terms_of_BSSN as AitoB
        AitoB.ADM_in_terms_of_BSSN()

        gammaDD = AitoB.gammaDD
        gammaDD_dD = AitoB.gammaDDdD
        KDD = AitoB.KDD
        
    elif formalism == "ADM":
        betaU = ixp.declarerank1("betaU") 
        betaU_dD = ixp.declarerank2("betaU_dD", "nosym")
        gammaDD = ixp.declarerank2("gammaDD","sym01")
        gammaDD_dD = ixp.declarerank3("gammaDD_dD","sym01")
        KDD = ixp.declarerank2("KDD","sym01")
        
    else:
        print("ERROR: "+formalism+" formalism not supported")
        
    # First define hydrodynamical quantities
    u4U = ixp.declarerank1("u4U", DIM=4)
    smallb4U = ixp.declarerank1("smallb4U", DIM=4)
    smallbsquared = sp.symbols("smallbsquared")
    ValenciavU = ixp.declarerank1("ValenciavU", DIM=3)
    VU = ixp.declarerank1("VU", DIM=3)
    rho_b, P, h, epsilon = sp.symbols('rhob P h epsilon',real=True)
    BU = ixp.declarerank1("BU", DIM=3)
    
    alpha = sp.symbols("alpha")
    alpha_dD = ixp.declarerank1("alpha_dD", DIM=3)

    

# Step 2.b.i: Define b^mu.
def compute_smallb4U(gammaDD,betaU,alpha, u4U, BU, sqrt4pi):
    global smallb4U
    import BSSN.ADMBSSN_tofrom_4metric as AB4m
    AB4m.g4DD_ito_BSSN_or_ADM("ADM",gammaDD,betaU,alpha)

    u4D = ixp.zerorank1(DIM=4)
    for mu in range(4):
        for nu in range(4):
            u4D[mu] += AB4m.g4DD[mu][nu]*u4U[nu]
    smallb4U = ixp.zerorank1(DIM=4)
    u4_dot_B_notilde = sp.sympify(0)
    for i in range(3):
        u4_dot_B_notilde += sp.simplify(u4D[i+1]*BU[i]).doit()
#         print(fr"finished {i} of simplify")

    # b^0 = (u_j B^j)/[alpha * sqrt(4 pi)]
    smallb4U[0] = sp.simplify(u4_dot_B_notilde / (alpha*sqrt4pi)).doit()
#     print(fr"finished smallb4U0 simplify")
    # b^i = [B^i + (u_j B^j) u^i]/[alpha * u^0 * sqrt(4 pi)]
    for i in range(3):
        smallb4U[i+1] = sp.simplify((BU[i] + u4_dot_B_notilde*u4U[i+1]) / (alpha*u4U[0]*sqrt4pi)).doit()


# Step 2.c: Define b^2.
def compute_smallbsquared(gammaDD,betaU,alpha, smallb4U):
    global smallbsquared
    import BSSN.ADMBSSN_tofrom_4metric as AB4m
    AB4m.g4DD_ito_BSSN_or_ADM("ADM",gammaDD,betaU,alpha)

    smallbsquared = sp.sympify(0)
    for mu in range(4):
        for nu in range(4):
            smallbsquared += AB4m.g4DD[mu][nu]*smallb4U[mu]*smallb4U[nu]


# Step 2.b: Define T^{mu nu} (a 4-dimensional tensor)
def compute_T4UU(gammaDD,betaU,alpha, rho_b, P, h, u4U, smallb4U, smallbsquared):
    global T4UU

    # Then define g^{mu nu} in terms of the ADM quantities:
    import BSSN.ADMBSSN_tofrom_4metric as AB4m
    AB4m.g4UU_ito_BSSN_or_ADM("ADM",gammaDD,betaU,alpha)

    # Finally compute T^{mu nu}
    T4UU = ixp.zerorank2(DIM=4)
    for mu in range(4):
        for nu in range(4):
            T4UU[mu][nu] = (rho_b * h + smallbsquared) * u4U[mu]*u4U[nu]\
                         + (P + sp.Rational(1,2)*smallbsquared)*AB4m.g4UU[mu][nu]\
                         - smallb4U[mu]*smallb4U[nu]

# Step 3.b: Define T^{mu}_{nu} (a 4-dimensional tensor)
def compute_T4UD(gammaDD,betaU,alpha, T4UU):
    global T4UD
    # Next compute T^mu_nu = T^{mu delta} g_{delta nu}, needed for S_tilde flux.
    # First we'll need g_{alpha nu} in terms of ADM quantities:
    import BSSN.ADMBSSN_tofrom_4metric as AB4m
    AB4m.g4DD_ito_BSSN_or_ADM("ADM",gammaDD,betaU,alpha)
    T4UD = ixp.zerorank2(DIM=4)
    for mu in range(4):
        for nu in range(4):
            for delta in range(4):
                T4UD[mu][nu] += T4UU[mu][delta]*AB4m.g4DD[delta][nu]


# Step 3: Writing the conservative variables in terms of the primitive variables
def compute_sqrtgammaDET(gammaDD):
    global sqrtgammaDET
    _gammaUU, gammaDET = ixp.symm_matrix_inverter3x3(gammaDD) # _gammaUU unused.
    sqrtgammaDET = sp.sqrt(gammaDET)

def compute_rho_star(alpha, sqrtgammaDET, rho_b,u4U):
    global rho_star
    # Compute rho_star:
    rho_star = alpha*sqrtgammaDET*rho_b*u4U[0]

def compute_tau_tilde(alpha, sqrtgammaDET, T4UU,rho_star):
    global tau_tilde
    tau_tilde = alpha**2*sqrtgammaDET*T4UU[0][0] - rho_star

def compute_S_tildeD(alpha, sqrtgammaDET, T4UD):
    global S_tildeD
    S_tildeD = ixp.zerorank1(DIM=3)
    for i in range(3):
        S_tildeD[i] = alpha*sqrtgammaDET*T4UD[0][i+1]


# Step 4: Define the fluxes for the GRHD equations
# Step 4.a: vU from u4U may be needed for computing rho_star_flux from u4U
def compute_vU_from_u4U__no_speed_limit(u4U):
    global VU
    # Now compute v^i = u^i/u^0:
    VU = ixp.zerorank1(DIM=3)
    for j in range(3):
        VU[j] = u4U[j+1]/u4U[0]

# Step 4.b: rho_star flux
def compute_rho_star_fluxU(vU, rho_star):
    global rho_star_fluxU
    rho_star_fluxU = ixp.zerorank1(DIM=3)
    for j in range(3):
        rho_star_fluxU[j] = rho_star*vU[j]

# Step 4.c: tau_tilde flux
def compute_tau_tilde_fluxU(alpha, sqrtgammaDET, vU,T4UU,rho_star):
    global tau_tilde_fluxU
    tau_tilde_fluxU = ixp.zerorank1(DIM=3)
    for j in range(3):
        tau_tilde_fluxU[j] = alpha**2*sqrtgammaDET*T4UU[0][j+1] - rho_star*vU[j]

# Step 4.d: S_tilde flux
def compute_S_tilde_fluxUD(alpha, sqrtgammaDET, T4UD):
    global S_tilde_fluxUD
    S_tilde_fluxUD = ixp.zerorank2(DIM=3)
    for j in range(3):
        for i in range(3):
            S_tilde_fluxUD[j][i] = alpha*sqrtgammaDET*T4UD[j+1][i+1]

# Step 5: Define source terms on RHSs of GRHD equations
# Step 5.a: tau_tilde RHS source term s
def compute_tau_tilde_source_term(KDD,betaU,alpha, sqrtgammaDET,alpha_dD, T4UU):
    global tau_tilde_source_term
    tau_tilde_source_term = sp.sympify(0)
    for i in range(3):
        for j in range(3):
            tau_tilde_source_term += (T4UU[0][0]*betaU[i]*betaU[j] + 2*T4UU[0][i+1]*betaU[j] + T4UU[i+1][j+1])*KDD[i][j]

    for i in range(3):
        tau_tilde_source_term += -(T4UU[0][0]*betaU[i] + T4UU[0][i+1])*alpha_dD[i]

    tau_tilde_source_term *= alpha*sqrtgammaDET

# Step 5.b: Define source term on RHS of $\tilde{S}_i$ equation
# Step 5.b.i: Compute g_{mu nu, i}, needed for the S tilde source term
def compute_g4DD_zerotimederiv_dD(gammaDD,betaU,alpha, gammaDD_dD,betaU_dD,alpha_dD):
    global g4DD_zerotimederiv_dD
    # Eq. 2.121 in B&S
    betaD = ixp.zerorank1()
    for i in range(3):
        for j in range(3):
            betaD[i] += gammaDD[i][j]*betaU[j]

    # gammaDD_dD = ixp.declarerank3("gammaDD_dDD","sym12",DIM=3)
    # betaU_dD   = ixp.declarerank2("betaU_dD"   ,"nosym",DIM=3)
    betaDdD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            for k in range(3):
                # Recall that betaD[i] = gammaDD[i][j]*betaU[j] (Eq. 2.121 in B&S)
                betaDdD[i][k] += gammaDD_dD[i][j][k]*betaU[j] + gammaDD[i][j]*betaU_dD[j][k]

    # Eq. 2.122 in B&S
    g4DD_zerotimederiv_dD = ixp.zerorank3(DIM=4)
    for k in range(3):
        # Recall that g4DD[0][0] = -alpha^2 + betaU[j]*betaD[j]
        g4DD_zerotimederiv_dD[0][0][k+1] += -2*alpha*alpha_dD[k]
        for j in range(3):
            g4DD_zerotimederiv_dD[0][0][k+1] += betaU_dD[j][k]*betaD[j] + betaU[j]*betaDdD[j][k]

    for i in range(3):
        for k in range(3):
            # Recall that g4DD[i][0] = g4DD[0][i] = betaD[i]
            g4DD_zerotimederiv_dD[i+1][0][k+1] = g4DD_zerotimederiv_dD[0][i+1][k+1] = betaDdD[i][k]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                # Recall that g4DD[i][j] = gammaDD[i][j]
                g4DD_zerotimederiv_dD[i+1][j+1][k+1] = gammaDD_dD[i][j][k]

# Step 5.b.ii: S_tilde source terms
def compute_S_tilde_source_termD(alpha, sqrtgammaDET,g4DD_zerotimederiv_dD, T4UU):
    global S_tilde_source_termD
    S_tilde_source_termD = ixp.zerorank1(DIM=3)
    for i in range(3):
        for mu in range(4):
            for nu in range(4):
                S_tilde_source_termD[i] += sp.Rational(1,2)*alpha*sqrtgammaDET*T4UU[mu][nu]*g4DD_zerotimederiv_dD[mu][nu][i+1]

                
# Step 6.a: Convert Valencia 3-velocity v_{(n)}^i into u^\mu, and apply a speed limiter
#           Speed-limited ValenciavU is output to rescaledValenciavU global.
def u4U_in_terms_of_ValenciavU__rescale_ValenciavU_by_applying_speed_limit(alpha,betaU,gammaDD, ValenciavU):
    # Inputs:  Metric lapse alpha, shift betaU, 3-metric gammaDD, Valencia 3-velocity ValenciavU
    # Outputs (as globals): u4U_ito_ValenciavU, rescaledValenciavU

    # R = gamma_{ij} v^i v^j
    R = sp.sympify(0)
    for i in range(3):
        for j in range(3):
            R += gammaDD[i][j]*ValenciavU[i]*ValenciavU[j]
    
    Rmax = 1 - 1 / (GAMMA_SPEED_LIMIT * GAMMA_SPEED_LIMIT)
    # Now, we set Rstar = min(Rmax,R):
    # If R <  Rmax, then Rstar = 0.5*(Rmax+R-Rmax+R) = R
    # If R >= Rmax, then Rstar = 0.5*(Rmax+R+Rmax-R) = Rmax
    Rstar = sp.Rational(1, 2) * (Rmax + R - nrpyAbs(Rmax - R))

    # We add TINYDOUBLE to R below to avoid a 0/0, which occurs when
    #    ValenciavU == 0 for all Valencia 3-velocity components.
    # "Those tiny *doubles* make me warm all over
    #  with a feeling that I'm gonna love you till the end of time."
    #    - Adapted from Connie Francis' "Tiny Bubbles"
    TINYDOUBLE = par.Cparameters("#define",thismodule,"TINYDOUBLE",1e-100)

    # The rescaled (speed-limited) Valencia 3-velocity
    #     is given by, v_{(n)}^i = sqrt{Rstar/R} v^i
    global rescaledValenciavU
    rescaledValenciavU = ixp.zerorank1(DIM=3)
    for i in range(3):
        # If R == 0, then Rstar == 0, so sqrt( Rstar/(R+TINYDOUBLE) )=sqrt(0/1e-100) = 0
        #   If your velocities are of order 1e-100 and this is physically
        #   meaningful, there must be something wrong with your unit conversion.
        rescaledValenciavU[i] = ValenciavU[i]*sp.sqrt(Rstar/(R + TINYDOUBLE))

    # Finally compute u^mu in terms of Valenciav^i
    # u^0 = 1/(alpha-sqrt(1-R^*))
    global u4U_ito_ValenciavU
    u4U_ito_ValenciavU = ixp.zerorank1(DIM=4)
    u4U_ito_ValenciavU[0] = 1/(alpha*sp.sqrt(1-Rstar))
    # u^i = u^0 ( alpha v^i_{(n)} - beta^i ), where v^i_{(n)} is the Valencia 3-velocity
    for i in range(3):
        u4U_ito_ValenciavU[i+1] = u4U_ito_ValenciavU[0] * (alpha * rescaledValenciavU[i] - betaU[i])

# Step 6.b: Convert v^i into u^\mu, and apply a speed limiter.
#           Speed-limited vU is output to rescaledvU global.
def u4U_in_terms_of_vU__rescale_vU_by_applying_speed_limit(alpha,betaU,gammaDD, VU):
    ValenciavU = ixp.zerorank1(DIM=3)
    for i in range(3):
        ValenciavU[i] = (VU[i] + betaU[i])/alpha
    u4U_in_terms_of_ValenciavU__rescale_ValenciavU_by_applying_speed_limit(alpha,betaU,gammaDD, ValenciavU)
    # Since ValenciavU is written in terms of vU,
    #   u4U_ito_ValenciavU is actually u4U_ito_vU
    global u4U_ito_VU
    u4U_ito_VU = ixp.zerorank1(DIM=4)
    for mu in range(4):
        u4U_ito_VU[mu] = u4U_ito_ValenciavU[mu]
    # Finally compute the rescaled (speed-limited) vU
    global rescaledVU
    rescaledVU = ixp.zerorank1(DIM=3)
    for i in range(3):
        rescaledVU[i] = alpha * rescaledValenciavU[i] - betaU[i]
    
