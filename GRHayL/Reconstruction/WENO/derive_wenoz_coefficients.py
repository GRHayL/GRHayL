#!/usr/bin/env python3
"""
Generate the finite-volume WENO-Z tables used in the higher-order GRHayL
reconstruction routines.

Conventions
-----------
- The cell width is normalized to 1.
- The center cell is at x = 0.
- The interface to reconstruct is the right face x = 1/2.
- The r-th order candidate substencil S_k contains the relative cells
  {k-(r-1), ..., k}.

For each substencil, this script derives:
- alpha_{k,j}: the coefficient multiplying the j-th cell average on substencil
  k when evaluating the reconstruction polynomial at x = 1/2.
- gamma_k: the optimal linear weights that reproduce the full degree-(2r-2)
  polynomial built from the complete 2r-1 cell stencil.
- B_k: the Jiang-Shu smoothness matrix such that beta_k = q_k^T B_k q_k.

The higher-order WENO-Z tau coefficients are the standard Castro-Costa-Don
choices:
  r = 4: [1, 3, -3, -1]
  r = 5: [1, 2, -6, 2, 1]
They are chosen so that the smooth-region Taylor expansion of sum_k c_k beta_k
has higher-order cancellation.
"""

import sympy as sp

x = sp.symbols("x")

TAU_COEFFS = {
    4: [1, 3, -3, -1],
    5: [1, 2, -6, 2, 1],
}


def cell_interval(j):
    return (sp.Rational(2 * j - 1, 2), sp.Rational(2 * j + 1, 2))


def relative_offsets(r, k):
    return [k - (r - 1) + j for j in range(r)]


def candidate_polynomial(r, k):
    q = sp.symbols(f"q0:{r}")
    a = sp.symbols(f"a0:{r}")
    p = sum(a[m] * x**m for m in range(r))
    eqs = []

    for idx, offset in enumerate(relative_offsets(r, k)):
        xl, xr = cell_interval(offset)
        eqs.append(sp.Eq(sp.integrate(p, (x, xl, xr)), q[idx]))

    sol = sp.solve(eqs, a, dict=True)[0]
    return sp.expand(p.subs(sol)), q


def candidate_matrix(r):
    alpha = []

    for k in range(r):
        poly, q = candidate_polynomial(r, k)
        face_value = sp.expand(poly.subs(x, sp.Rational(1, 2)))
        alpha.append([sp.simplify(face_value.coeff(qj)) for qj in q])

    return sp.Matrix(alpha)


def full_stencil_face_coeffs(r):
    q = sp.symbols(f"q0:{2 * r - 1}")
    a = sp.symbols(f"a0:{2 * r - 1}")
    p = sum(a[m] * x**m for m in range(2 * r - 1))
    eqs = []

    for idx, offset in enumerate(range(-(r - 1), r)):
        xl, xr = cell_interval(offset)
        eqs.append(sp.Eq(sp.integrate(p, (x, xl, xr)), q[idx]))

    sol = sp.solve(eqs, a, dict=True)[0]
    face_value = sp.expand(p.subs(sol).subs(x, sp.Rational(1, 2)))
    return [sp.simplify(face_value.coeff(qj)) for qj in q]


def optimal_weights(r, alpha=None):
    if alpha is None:
        alpha = candidate_matrix(r)

    full_coeffs = full_stencil_face_coeffs(r)
    gamma = sp.symbols(f"g0:{r}")
    eqs = []

    for n in range(2 * r - 1):
        coeff = 0
        for k in range(r):
            local_idx = n - k
            if 0 <= local_idx < r:
                coeff += gamma[k] * alpha[k, local_idx]
        eqs.append(sp.Eq(coeff, full_coeffs[n]))

    sol = list(sp.linsolve(eqs, gamma))[0]
    return sp.Matrix(sol)


def smoothness_matrices(r):
    beta_mats = []

    for k in range(r):
        poly, q = candidate_polynomial(r, k)
        beta = 0

        for m in range(1, r):
            beta += sp.integrate(
                sp.diff(poly, x, m) ** 2,
                (x, sp.Rational(-1, 2), sp.Rational(1, 2)),
            )

        beta_mats.append(sp.simplify(sp.hessian(sp.expand(beta), q) / 2))

    return beta_mats


def to_c(expr):
    return sp.ccode(sp.simplify(expr))


def print_table(name, data):
    if isinstance(data, sp.Matrix):
        rows = data.tolist()
    else:
        rows = data

    print(f"{name} = {{")
    for row in rows:
        if isinstance(row, (list, tuple)):
            print("  {" + ", ".join(to_c(val) for val in row) + "},")
        else:
            print("  " + to_c(row) + ",")
    print("}")


def emit(r):
    alpha = candidate_matrix(r)
    gamma = optimal_weights(r, alpha)
    beta = smoothness_matrices(r)

    print(f"r = {r}")
    print_table("alpha", alpha)
    print_table("gamma", [list(gamma)])
    print_table("tau", [TAU_COEFFS[r]])
    print("beta = {")
    for mat in beta:
      print("  {")
      for row in mat.tolist():
        print("    {" + ", ".join(to_c(val) for val in row) + "},")
      print("  },")
    print("}")


if __name__ == "__main__":
    emit(4)
    print()
    emit(5)
