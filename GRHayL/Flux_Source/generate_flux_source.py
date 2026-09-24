"""
Generate GRHayL GRMHD source terms, characteristic speeds, and HLLE fluxes.

The source terms use the ADM form of the GRMHD evolution equations. Face
fluxes construct the Cartesian GRMHD tensor from NRPy 2 metric and magnetic
expressions and reconstructed face states. Run through generate_flux_source.sh
to use the pinned NRPy 2 commit.

Authors: GRHayL contributors
"""

import argparse
import os
import re
from pathlib import Path
from typing import List, NamedTuple, Set, Tuple

import nrpy
import nrpy.c_codegen as ccg
import nrpy.indexedexp as ixp
import sympy as sp
from nrpy.equations.general_relativity.g4munu_conversions import (
    ADM_to_g4DD,
    ADM_to_g4UU,
)
from nrpy.equations.grhd.characteristic_speeds import (
    find_cmax_cmin as find_grhd_cmax_cmin,
)
from nrpy.equations.grmhd.characteristic_speeds import compute_v02
from nrpy.equations.grmhd.GRMHD_equations import compute_smallb2

# Step P1: Confirm the shell script selected the pinned NRPy 2 checkout.
nrpy_root = os.environ.get("NRPY_ROOT")
if not nrpy_root or not (Path(nrpy_root) / "nrpy" / "__init__.py").is_file():
    raise RuntimeError("NRPY_ROOT must name the checked-out NRPy 2 repository")
if (
    Path(nrpy.__file__).resolve()
    != (Path(nrpy_root) / "nrpy" / "__init__.py").resolve()
):
    raise RuntimeError("Python imported NRPy outside NRPY_ROOT")

SOURCE_DIR = Path(__file__).resolve().parent
VARIANTS = ("hybrid", "hybrid_entropy", "tabulated", "tabulated_entropy")


class FaceSymbols(NamedTuple):
    """Symbolic ADM metric and primitive variables at one cell face."""

    alpha: sp.Expr
    beta: List[sp.Expr]
    gamma: List[List[sp.Expr]]
    u_r: List[sp.Expr]
    u_l: List[sp.Expr]
    B_r: List[sp.Expr]
    B_l: List[sp.Expr]
    P_r: sp.Expr
    P_l: sp.Expr
    rho_r: sp.Expr
    rho_l: sp.Expr
    reads: List[str]


class FluxState(NamedTuple):
    """Conserved variables U and physical fluxes F for one face state."""

    rho_cons: sp.Expr
    rho_flux: sp.Expr
    ye_cons: sp.Expr
    ye_flux: sp.Expr
    entropy_cons: sp.Expr
    entropy_flux: sp.Expr
    tau_cons: sp.Expr
    tau_flux: sp.Expr
    momentum_cons: List[sp.Expr]
    momentum_flux: List[sp.Expr]


def code(expressions: List[sp.Expr], outputs: List[str]) -> str:
    """
    Convert GRMHD expressions to C with NRPy 2 common-subexpression elimination.

    Golden Kernels replaces small integer powers with products or reciprocals.

    :param expressions: Symbolic source terms, speeds, or fluxes.
    :param outputs: C lvalues receiving the corresponding expressions.
    :return: C declarations and assignments for the equations.
    """
    return ccg.c_codegen(
        expressions,
        outputs,
        fp_type="double",
        fp_type_alias="double",
        verbose=False,
        include_braces=False,
        enable_GoldenKernels=True,
        enable_cse_preprocess=True,
        cse_sorting="canonical",
    )


def symbols(prefix: str, length: int) -> List[sp.Expr]:
    """Declare components of a rank-one NRPy indexed expression.

    :param prefix: Base name of the component symbols.
    :param length: Number of components.
    :return: Component symbols in index order.
    """
    return ixp.declarerank1(prefix, dimension=length)


def metric_symbols(
    suffix: str = "",
) -> Tuple[sp.Expr, List[sp.Expr], List[List[sp.Expr]]]:
    """Declare ADM lapse, shift, and symmetric spatial metric components.

    :param suffix: Name suffix distinguishing a face metric from a cell metric.
    :return: Lapse, shift, and symmetric spatial metric components.
    """
    alpha = sp.Symbol("alpha" + suffix)
    beta = symbols("beta" + suffix + "U", 3)
    gamma = ixp.declarerank2("gamma" + suffix + "DD", symmetry="sym01", dimension=3)
    return alpha, beta, gamma


def metric_reads(
    pointer: str, alpha: sp.Expr, beta: List[sp.Expr], gamma: List[List[sp.Expr]]
) -> List[str]:
    """Read one GRHayL ADM metric into NRPy's scalar component names.

    :param pointer: C pointer to the metric structure.
    :param alpha: Symbol receiving the lapse.
    :param beta: Symbols receiving the shift components.
    :param gamma: Symbols receiving the spatial metric components.
    :return: C declarations reading lapse, shift, and unique metric entries.
    """
    lines = [f"const double {alpha} = {pointer}->lapse;"]
    lines += [f"const double {beta[i]} = {pointer}->betaU[{i}];" for i in range(3)]
    lines += [
        f"const double {gamma[i][j]} = {pointer}->gammaDD[{i}][{j}];"
        for i in range(3)
        for j in range(i, 3)
    ]
    return lines


def primitive_reads(
    pointer: str, side: str = ""
) -> Tuple[List[sp.Expr], List[sp.Expr], sp.Expr, sp.Expr, List[str]]:
    """
    Read a primitive state and form u^i = v^i u^0.

    GRHayL stores magnetic components after division by sqrt(4 pi), matching
    NRPy 2's BmagU convention. The caller has already limited the velocity.

    :param pointer: C pointer to a GRHayL primitive state.
    :param side: Face suffix, either r or l; empty for cell-centered source terms.
    :return: Four-velocity, magnetic field, pressure, density, and C reads.
    """
    u = symbols("u4" + side + "U" if side else "u4U", 4)
    B = symbols("B" + side + "U" if side else "BU", 3)
    pressure = sp.Symbol("P_" + side if side else "P")
    density = sp.Symbol("rhob_" + side if side else "rhob")
    lines = [f"const double {u[0]} = {pointer}->u0;"]
    lines += [f"const double {u[i + 1]} = {pointer}->vU[{i}]*{u[0]};" for i in range(3)]
    lines += [f"const double {B[i]} = {pointer}->BU[{i}];" for i in range(3)]
    lines += [
        f"const double {pressure} = {pointer}->press;",
        f"const double {density} = {pointer}->rho;",
    ]
    return u, B, pressure, density, lines


def eos_reads(two_states: bool) -> str:
    """
    Call the GRHayL EOS before reading mutable primitive values.

    A tabulated EOS can change density, pressure, and electron fraction. A
    failed callback leaves conservative outputs untouched.

    :param two_states: Whether the C function receives right and left states.
    :return: C declarations and checked EOS calls.
    """
    if two_states:
        return """  // The EOS may update both face states before the NRPy equations read them.
  double h_r, h_l, cs2_r, cs2_l;
  ghl_error_codes_t error = ghl_compute_h_and_cs2(eos, prims_r, &h_r, &cs2_r);
  if(error != ghl_success) return error;
  error = ghl_compute_h_and_cs2(eos, prims_l, &h_l, &cs2_l);
  if(error != ghl_success) return error;
"""
    return """  // The EOS may update the cell state before the source equations read it.
  double h, cs2;
  const ghl_error_codes_t error = ghl_compute_h_and_cs2(eos, prims, &h, &cs2);
  if(error != ghl_success) return error;
"""


def write_function(
    destination: Path,
    name: str,
    params: str,
    preamble: str,
    body: str,
    arguments: str,
    equation_comment: str,
) -> None:
    """
    Write checked C function and legacy aborting wrapper.

    :param destination: Generated C source file.
    :param name: Public legacy function name.
    :param params: C parameter declarations shared by both functions.
    :param preamble: C input reads and validation before equation assignments.
    :param body: NRPy-generated equation assignments.
    :param arguments: C arguments passed from the legacy wrapper.
    :param equation_comment: Physics equation evaluated by the checked C function.
    """
    checked = name + "_checked"
    # Step 1: NRPy 2 emits equation declarations at column zero. Indent those
    # lines in the C function without changing nested wave-speed checks.
    equations = "".join(
        "  " + line if line.strip() and not line[0].isspace() else line
        for line in (preamble + body).splitlines(keepends=True)
    )
    # Step 2: Keep checked error handling and the legacy aborting wrapper paired.
    content = (
        "// Generated by generate_flux_source.py using the pinned NRPy 2 checkout.\n"
        '#include "ghl_flux_source.h"\n\n'
        f"ghl_error_codes_t {checked}({params}) {{\n"
        f"  // {equation_comment}\n"
        f"{equations}  return ghl_success;\n"
        "}\n\n"
        f"void {name}({params}) {{\n"
        f"  ghl_abort_if_error({checked}({arguments}));\n"
        "}\n"
    )
    destination.write_text(content)


def grmhd_stress_energy(
    alpha: sp.Expr,
    beta: List[sp.Expr],
    gamma: List[List[sp.Expr]],
    u: List[sp.Expr],
    B: List[sp.Expr],
    density: sp.Expr,
    pressure: sp.Expr,
    enthalpy: sp.Expr,
) -> List[List[sp.Expr]]:
    """
    Form the total GRMHD stress-energy tensor in the ADM coordinate basis.

    NRPy 2 supplies the ADM metric and comoving magnetic norm. Combining
    rho*h with b^2 before multiplication by u^mu*u^nu preserves precision
    in momentum fluxes where magnetic and fluid terms nearly cancel.

    :param alpha: ADM lapse.
    :param beta: ADM shift components.
    :param gamma: Covariant spatial metric components.
    :param u: Contravariant fluid four-velocity.
    :param B: Eulerian magnetic field divided by sqrt(4 pi).
    :param density: Baryon density.
    :param pressure: Fluid pressure.
    :param enthalpy: Specific enthalpy.
    :return: Contravariant total stress-energy tensor.
    """
    # Step 1: Contract u_i B^i, grouped by B^i for CSE and C evaluation.
    g4DD = ADM_to_g4DD(gamma, beta, alpha)
    u_dot_B = sum(
        B[i] * sum(g4DD[i + 1][mu] * u[mu] for mu in range(4)) for i in range(3)
    )
    # Step 2: Form b^mu and b^2 using the NRPy 2 GRMHD contraction.
    smallb4U = [u_dot_B / alpha] + [
        (B[i] + u_dot_B * u[i + 1]) / (alpha * u[0]) for i in range(3)
    ]
    smallb2 = compute_smallb2(gamma, beta, alpha, smallb4U)
    g4UU = ADM_to_g4UU(gamma, beta, alpha)
    # Step 3: Assemble the fluid and magnetic terms in the Valencia tensor.
    return [
        [
            (density * enthalpy + smallb2) * u[mu] * u[nu]
            + (pressure + sp.Rational(1, 2) * smallb2) * g4UU[mu][nu]
            - smallb4U[mu] * smallb4U[nu]
            for nu in range(4)
        ]
        for mu in range(4)
    ]


def source_expressions() -> Tuple[List[sp.Expr], List[str]]:
    """
    Build ADM source terms for S_i and tau from the GRMHD stress tensor.

    The caller supplies spatial metric derivatives and extrinsic curvature;
    no finite-difference stencil is generated here. The formula for S_i is
    alpha sqrt(gamma) T^{mu nu} partial_i g_{mu nu}/2.

    :return: S_0, S_1, S_2, tau expressions and matching C input reads.
    """
    alpha, beta, gamma = metric_symbols()
    u, B, pressure, density, reads = primitive_reads("prims")
    K = ixp.declarerank2("KDD", symmetry="sym01", dimension=3)
    alpha_dD = symbols("alpha_dD", 3)
    beta_dD = [[sp.Symbol(f"betaU_dD{i}{d}") for d in range(3)] for i in range(3)]
    gamma_dD = [
        [
            [sp.Symbol(f"gammaDD_dD{min(i,j)}{max(i,j)}{d}") for d in range(3)]
            for j in range(3)
        ]
        for i in range(3)
    ]
    volume = sp.sqrt(sp.det(sp.Matrix(gamma)))

    # Step 1: Use the same GRMHD tensor for source and face-flux generation.
    tensor = grmhd_stress_energy(
        alpha, beta, gamma, u, B, density, pressure, sp.Symbol("h")
    )

    # Step 2: Contract T^{mu nu} with K_ij and partial_i alpha for tau source.
    tau = sum(
        (
            tensor[0][0] * beta[i] * beta[j]
            + 2 * tensor[0][i + 1] * beta[j]
            + tensor[i + 1][j + 1]
        )
        * K[i][j]
        for i in range(3)
        for j in range(3)
    )
    tau -= sum(
        (tensor[0][0] * beta[i] + tensor[0][i + 1]) * alpha_dD[i] for i in range(3)
    )
    tau *= alpha * volume

    # Step 3: Form g_{mu nu,i} from ADM derivatives. Reuse beta_i and beta_{i,j}
    # in g_{00,i} and g_{0i,j}; S_i = alpha sqrt(gamma) T^{mu nu} g_{mu nu,i}/2.
    betaD = [sum(gamma[i][j] * beta[j] for j in range(3)) for i in range(3)]
    betaD_dD = [
        [
            sum(
                gamma_dD[i][j][d] * beta[j] + gamma[i][j] * beta_dD[j][d]
                for j in range(3)
            )
            for d in range(3)
        ]
        for i in range(3)
    ]
    momentum = []
    for d in range(3):
        g4DD_dD = [[sp.sympify(0) for _ in range(4)] for _ in range(4)]
        g4DD_dD[0][0] = -2 * alpha * alpha_dD[d] + sum(
            beta_dD[j][d] * betaD[j] + beta[j] * betaD_dD[j][d] for j in range(3)
        )
        for i in range(3):
            g4DD_dD[0][i + 1] = g4DD_dD[i + 1][0] = betaD_dD[i][d]
            for j in range(3):
                g4DD_dD[i + 1][j + 1] = gamma_dD[i][j][d]
        momentum.append(
            sp.Rational(1, 2)
            * alpha
            * volume
            * sum(
                tensor[mu][nu] * g4DD_dD[mu][nu] for mu in range(4) for nu in range(4)
            )
        )

    reads += metric_reads("metric", alpha, beta, gamma)
    reads += [
        f"const double {K[i][j]} = curv->K[{i}][{j}];"
        for i in range(3)
        for j in range(i, 3)
    ]
    for d, direction in enumerate("xyz"):
        pointer = f"metric_derivs_{direction}"
        reads.append(f"const double {alpha_dD[d]} = {pointer}->lapse;")
        reads += [
            f"const double {beta_dD[i][d]} = {pointer}->betaU[{i}];" for i in range(3)
        ]
        reads += [
            f"const double {gamma_dD[i][j][d]} = {pointer}->gammaDD[{i}][{j}];"
            for i in range(3)
            for j in range(i, 3)
        ]
    return momentum + [tau], reads


def generate_source(destination: Path) -> None:
    """Write ``ghl_calculate_source_terms.c`` from ADM source expressions.

    :param destination: Empty output directory.
    """
    expressions, reads = source_expressions()
    params = (
        "const ghl_eos_parameters *restrict eos, ghl_primitive_quantities *restrict prims, "
        "const ghl_metric_quantities *restrict metric, "
        "const ghl_metric_quantities *restrict metric_derivs_x, "
        "const ghl_metric_quantities *restrict metric_derivs_y, "
        "const ghl_metric_quantities *restrict metric_derivs_z, "
        "const ghl_extrinsic_curvature *restrict curv, "
        "ghl_conservative_quantities *restrict cons"
    )
    arguments = "eos, prims, metric, metric_derivs_x, metric_derivs_y, metric_derivs_z, curv, cons"
    preamble = eos_reads(False) + "\n".join(reads) + "\n"
    body = code(expressions, [f"cons->SD[{i}]" for i in range(3)] + ["cons->tau"])
    write_function(
        destination / "ghl_calculate_source_terms.c",
        "ghl_calculate_source_terms",
        params,
        preamble,
        body,
        arguments,
        "Contract the GRMHD stress tensor with ADM metric derivatives and K_ij.",
    )


def face_symbols() -> FaceSymbols:
    """Declare right and left face primitives and their common ADM metric.

    :return: Symbolic face variables and C declarations reading their values.
    """
    alpha, beta, gamma = metric_symbols("_face")
    u_r, B_r, pressure_r, density_r, reads_r = primitive_reads("prims_r", "r")
    u_l, B_l, pressure_l, density_l, reads_l = primitive_reads("prims_l", "l")
    return FaceSymbols(
        alpha,
        beta,
        gamma,
        u_r,
        u_l,
        B_r,
        B_l,
        pressure_r,
        pressure_l,
        density_r,
        density_l,
        reads_r + reads_l + metric_reads("metric_face", alpha, beta, gamma),
    )


def wave_guard(direction: int) -> str:
    """Generate nonnegative, finite HLLE speed checks before EOS calls.

    :param direction: Coordinate direction of the face normal.
    :return: C checks and normalized HLLE wave-speed coefficients.
    """
    cmin, cmax = f"cmin_dirn{direction}", f"cmax_dirn{direction}"
    return f"""  // Validate HLLE wave speeds before calling the EOS or writing fluxes.
  const double wavespeed_scale = fmax(1.0, fmax(fabs({cmin}), fabs({cmax})));
  if(!isfinite({cmin}) || !isfinite({cmax})
     || {cmin} < -DBL_EPSILON * wavespeed_scale
     || {cmax} < -DBL_EPSILON * wavespeed_scale) return ghl_error_invalid_hlle_wavespeeds;
  const double cmin_floored = fmax({cmin}, 0.0);
  const double cmax_floored = fmax({cmax}, 0.0);
  if(cmin_floored > DBL_MAX - cmax_floored
     || (cmin_floored > 1.0 && cmax_floored > DBL_MAX / cmin_floored))
    return ghl_error_invalid_hlle_wavespeeds;
  const double wavespeed_sum = cmin_floored + cmax_floored;
  if(wavespeed_sum <= 0.0 || wavespeed_sum < 1.0 / DBL_MAX)
    return ghl_error_invalid_hlle_wavespeeds;
  const double cmin_weight = cmin_floored / wavespeed_sum;
  const double cmax_weight = cmax_floored / wavespeed_sum;
  const double dissipation_speed = cmin_floored * cmax_floored / wavespeed_sum;
"""


def generate_speeds(destination: Path, face: FaceSymbols) -> None:
    """
    Write characteristic speeds for all three face directions.

    NRPy 2 computes v_0^2 from b^2 and sound speed. Its GRHD root calculation
    then accepts v_0^2 for the GRMHD fast-wave estimate.

    :param destination: Empty output directory.
    :param face: Symbolic face metric and primitive states.
    """
    alpha, beta, gamma = face.alpha, face.beta, face.gamma
    u_r, u_l, B_r, B_l = face.u_r, face.u_l, face.B_r, face.B_l
    rho_r, rho_l, reads = face.rho_r, face.rho_l, face.reads
    h_r, h_l, cs2_r, cs2_l = sp.symbols("h_r h_l cs2_r cs2_l")
    v02_r = compute_v02(gamma, beta, alpha, u_r, B_r, rho_r, h_r, cs2_r)
    v02_l = compute_v02(gamma, beta, alpha, u_l, B_l, rho_l, h_l, cs2_l)
    # Step 1: Generate v_0^2 separately for the two face states. This prevents
    # CSE from expanding the magnetic contractions inside all four roots.
    speed_body = (
        "  // Magnetic fast-wave speed squared for right and left states.\n"
        "  double v02_r, v02_l;\n  {\n"
        + code([v02_r, v02_l], ["v02_r", "v02_l"])
        + "  }\n"
    )
    params = (
        "ghl_primitive_quantities *restrict prims_r, ghl_primitive_quantities *restrict prims_l, "
        "const ghl_eos_parameters *restrict eos, const ghl_metric_quantities *restrict metric_face, "
    )
    # Step 2: Solve the characteristic quadratic and write cmin, cmax.
    for direction in range(3):
        cmin, cmax = find_grhd_cmax_cmin(
            direction,
            gamma,
            beta,
            alpha,
            u_r,
            u_l,
            sp.Symbol("v02_r"),
            sp.Symbol("v02_l"),
        )
        name = f"ghl_calculate_characteristic_speed_dirn{direction}"
        full_params = (
            params + f"double *cmin_dirn{direction}, double *cmax_dirn{direction}"
        )
        arguments = f"prims_r, prims_l, eos, metric_face, cmin_dirn{direction}, cmax_dirn{direction}"
        speed_reads = [line for line in reads if not line.startswith("const double P_")]
        preamble = eos_reads(True) + "\n".join(speed_reads) + "\n"
        body = speed_body + code(
            [cmin, cmax], [f"*cmin_dirn{direction}", f"*cmax_dirn{direction}"]
        )
        write_function(
            destination / f"{name}.c",
            name,
            full_params,
            preamble,
            body,
            arguments,
            "Solve the GRMHD fast-wave characteristic quadratic on this face.",
        )


def state_fluxes(direction: int, face: FaceSymbols) -> List[FluxState]:
    """
    Compute U and physical F for each side of one face.

    The total tensor includes NRPy 2 magnetic stress-energy. Lowering one
    index directly from T^{mu nu} avoids extra magnetic contractions and keeps
    generated C speed close to the original GRHayL flux functions.

    :param direction: Coordinate direction of the face normal.
    :param face: Symbolic face metric and primitive states.
    :return: Right and left conserved variables and physical fluxes.
    """
    alpha, beta, gamma = face.alpha, face.beta, face.gamma
    volume = sp.sqrt(sp.det(sp.Matrix(gamma)))
    result: List[FluxState] = []
    for side, u, B, density, pressure in (
        ("r", face.u_r, face.B_r, face.rho_r, face.P_r),
        ("l", face.u_l, face.B_l, face.rho_l, face.P_l),
    ):
        # Step 1: Use the source-term tensor for each reconstructed face state.
        tensor = grmhd_stress_energy(
            alpha, beta, gamma, u, B, density, pressure, sp.Symbol(f"h_{side}")
        )
        # Step 2: Lower the second index for momentum density and flux.
        g4DD = ADM_to_g4DD(gamma, beta, alpha)
        mixed = [
            [
                sum(tensor[mu][sigma] * g4DD[sigma][nu] for sigma in range(4))
                for nu in range(4)
            ]
            for mu in range(4)
        ]
        transport = u[direction + 1] / u[0]
        density_cons = alpha * volume * density * u[0]
        density_flux = density_cons * transport
        ye = sp.Symbol(f"Y_e_{side}")
        entropy = sp.Symbol(f"S_{side}")
        ye_cons = ye * density_cons
        entropy_cons = alpha * volume * entropy * u[0]
        energy_cons = alpha**2 * volume * tensor[0][0] - density_cons
        energy_flux = alpha**2 * volume * tensor[0][direction + 1] - density_flux
        momentum_cons = [alpha * volume * mixed[0][i + 1] for i in range(3)]
        momentum_flux = [alpha * volume * mixed[direction + 1][i + 1] for i in range(3)]
        # Step 3: Return Valencia conserved fields and their physical fluxes.
        result.append(
            FluxState(
                density_cons,
                density_flux,
                ye_cons,
                ye_cons * transport,
                entropy_cons,
                entropy_cons * transport,
                energy_cons,
                energy_flux,
                momentum_cons,
                momentum_flux,
            )
        )
    return result


def generate_fluxes(destination: Path, face: FaceSymbols) -> None:
    """
    Write twelve HLLE flux functions across three directions and four EOS modes.

    :param destination: Empty output directory with variant subdirectories.
    :param face: Symbolic face metric and primitive states.
    """
    reads = face.reads
    params = (
        "ghl_primitive_quantities *restrict prims_r, ghl_primitive_quantities *restrict prims_l, "
        "const ghl_eos_parameters *restrict eos, const ghl_metric_quantities *restrict metric_face, "
    )
    cmin_weight, cmax_weight, dissipation_speed = sp.symbols(
        "cmin_weight cmax_weight dissipation_speed"
    )

    def hlle(
        right_flux: sp.Expr,
        left_flux: sp.Expr,
        right_cons: sp.Expr,
        left_cons: sp.Expr,
    ) -> sp.Expr:
        """Combine conserved variables and fluxes with checked HLLE weights.

        :param right_flux: Physical flux from the right state.
        :param left_flux: Physical flux from the left state.
        :param right_cons: Conserved variable from the right state.
        :param left_cons: Conserved variable from the left state.
        :return: HLLE flux for one conserved field.
        """
        return (
            cmin_weight * right_flux
            + cmax_weight * left_flux
            - dissipation_speed * (right_cons - left_cons)
        )

    # Step 1: Calculate physical fluxes once per direction for all EOS modes.
    for direction in range(3):
        right, left = state_fluxes(direction, face)
        outputs = [
            hlle(
                right.momentum_flux[i],
                left.momentum_flux[i],
                right.momentum_cons[i],
                left.momentum_cons[i],
            )
            for i in range(3)
        ]
        outputs += [
            hlle(right.rho_flux, left.rho_flux, right.rho_cons, left.rho_cons),
            hlle(right.tau_flux, left.tau_flux, right.tau_cons, left.tau_cons),
            hlle(
                right.entropy_flux,
                left.entropy_flux,
                right.entropy_cons,
                left.entropy_cons,
            ),
            hlle(right.ye_flux, left.ye_flux, right.ye_cons, left.ye_cons),
        ]
        # Step 2: Include electron fraction and entropy only for evolved modes.
        for variant in VARIANTS:
            names = [f"cons->SD[{i}]" for i in range(3)] + ["cons->rho", "cons->tau"]
            expressions = outputs[:5]
            extra_reads = []
            if "entropy" in variant:
                names.append("cons->entropy")
                expressions.append(outputs[5])
                extra_reads += [
                    f"const double S_{side} = prims_{side}->entropy;"
                    for side in ("r", "l")
                ]
            if "tabulated" in variant:
                names.append("cons->Y_e")
                expressions.append(outputs[6])
                extra_reads += [
                    f"const double Y_e_{side} = prims_{side}->Y_e;"
                    for side in ("r", "l")
                ]
            name = f"ghl_calculate_HLLE_fluxes_dirn{direction}_{variant}"
            full_params = (
                params
                + f"const double cmin_dirn{direction}, const double cmax_dirn{direction}, ghl_conservative_quantities *restrict cons"
            )
            arguments = f"prims_r, prims_l, eos, metric_face, cmin_dirn{direction}, cmax_dirn{direction}, cons"
            preamble = (
                wave_guard(direction)
                + eos_reads(True)
                + "\n".join(reads + extra_reads)
                + "\n"
            )
            body = code(expressions, names)
            path = destination / variant / f"{name}.c"
            write_function(
                path,
                name,
                full_params,
                preamble,
                body,
                arguments,
                "Combine the left and right Valencia states with the HLLE flux.",
            )


def manifest_outputs() -> Set[Path]:
    """Read names of generated C files from the Flux_Source build manifests.

    :return: C paths relative to ``GRHayL/Flux_Source``.
    """
    paths: Set[Path] = set()
    for directory in (SOURCE_DIR, *(SOURCE_DIR / variant for variant in VARIANTS)):
        manifest = (directory / "make.code.defn").read_text()
        paths.update(
            (directory.relative_to(SOURCE_DIR) / name)
            for name in re.findall(r"\b[A-Za-z0-9_]+\.c\b", manifest)
        )
    return paths


def main() -> None:
    """Generate the complete C source set in a new or empty directory.

    :raises RuntimeError: If generated C names differ from the build manifests.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_directory", type=Path)
    args = parser.parse_args()
    destination = args.output_directory.expanduser().resolve()
    if destination == SOURCE_DIR or SOURCE_DIR in destination.parents:
        parser.error("output directory must be outside GRHayL/Flux_Source")
    if destination.exists() and (
        not destination.is_dir() or any(destination.iterdir())
    ):
        parser.error("output directory must be new or empty")
    destination.mkdir(parents=True, exist_ok=True)
    for variant in VARIANTS:
        (destination / variant).mkdir()
    face = face_symbols()
    generate_source(destination)
    generate_speeds(destination, face)
    generate_fluxes(destination, face)
    generated = {path.relative_to(destination) for path in destination.rglob("*.c")}
    expected = manifest_outputs()
    if generated != expected:
        raise RuntimeError(
            f"generated output differs from build manifests: missing={expected - generated}, extra={generated - expected}"
        )


if __name__ == "__main__":
    main()
