# CompOSE to GRHayL regularized table converter

`compose_to_grhayl.py` converts one fixed CompOSE-generated HDF5 table into
the StellarCollapse HDF5 schema already consumed by GRHayL.  The supported
production profile is `sro-sly4-sna-141-regularized-v1`: CompOSE table 141,
SRO(SLy4) single-nucleus approximation, on its official 391 by 163 by 66
grid.

This is a regularized surrogate EOS.  It deliberately changes node values to
make the legacy GRHayL interpolation and temperature inversions usable.  It is
not the unmodified official table 141, a free-energy reconstruction, or a
claim that every regularized field remains an equilibrium SRO derivative.

## Input

Generate the table with CompOSE program banner 2.17, first-order
interpolation, and HDF output mode 2.  The input directory must provide these
four regular, nonsymlink files; other directory entries are not consumed:

```text
eoscompose.h5
eos.parameters
eos.quantities
eos.thermo
```

The profile fixes the grid, axes, controls, quantity identifiers, particle
masses, density and chemical-potential reference masses, and the relation
`Y_e=Y_q`.  The converter rejects filesystem symlinks, extra or reordered
identifiers, unexpected HDF5 objects or aliases, soft/external HDF5 links,
filters, layouts, ranks, or data types.  It does not download, parse raw
CompOSE tables, resample, or overwrite files.

## Run

```sh
python3 -m venv compose-venv
compose-venv/bin/python -m pip install -r tools/compose/requirements.txt
compose-venv/bin/python tools/compose/compose_to_grhayl.py \
  --input-dir GENERATED_DIR \
  --profile sro-sly4-sna-141-regularized-v1 \
  --output TABLE.h5
```

Success exits 0.  Command-line and profile errors exit 2.  Input validation,
regularization, verification, or publication failures exit 1.  Output is
written to an exclusive temporary file in the destination directory, closed,
reopened and fully validated, then published without replacing an existing
path.  Only converter-owned temporary files are removed after failure.

The result contains all 19 double-precision StellarCollapse fields with shape
`(N_Ye,N_T,N_rho)`, density fastest, and `have_rel_cs2=1`.  GRHayL must load it
with runtime sound-speed cleaning disabled.  The diagnostic
`/grhayl_compose/manifest_json` records fixed profile facts, policy, and
distortion measures; it contains no source hashes, dates, file times, or local
paths.

## Regularization contract

For each fixed `(Y_e,rho)` ray, equal-weight, minimum-slope isotonic
regression acts on `log10(eps+shift)` and `log10(P)`.  Entropy is integrated
from the regularized energy using a logarithmic-temperature mean.  Positive
temperature derivatives use a monotone PCHIP-style stencil; density
derivatives use second-order nonuniform stencils.

The pressure minimum spacing uses the magnitude of the loader-transformed
`ln(P_code)` and an eightfold safety margin over the loader's `1e-10` inverse
tolerance.  Reopened-output validation rejects adjacent pressure gaps that
the inverse routine could mistake for an exact initial-temperature match.

Writing `e_j` for regularized energy in MeV per baryon and `S_j` for entropy
in units of Boltzmann's constant per baryon, the entropy update is

```text
S_0 = max(S_raw(T_0),0)
T_L = (T_j-T_{j-1})/ln(T_j/T_{j-1})
S_j = S_{j-1} + (e_j-e_{j-1})/T_L
```

The two pressure-derivative contributions are jointly projected so that

```text
B = rho*dpdrhoe + (P/rho)*dpderho
1e-14 <= B/[rho*(c^2+eps)+P] <= 1-1e-12
cs2 = c^2*B/[rho*(c^2+eps)+P]
gamma = B/P
```

Fractions are minimum-L2 projected onto nonnegative baryon and charge closure
at every node.  Absent heavy nuclei use the inert `Abar=1,Zbar=0,Xh=0`
sentinel and the explicit three-species projection

```text
Xa = clamp((1-Xn_raw-Xp_raw+2*Xa_raw)/3, 0, 2*min(Ye,1-Ye))
Xp = Ye-Xa/2
Xn = 1-Ye-Xa/2
Xh = 0
```

Chemical fields are canonicalized to
`muhat=mu_n-mu_p` and `munu=mu_e-muhat`.

Independent trilinear interpolation of `Abar`, `Zbar`, and `Xh` does not
preserve their nonlinear charge product away from nodes.  The manifest reports
a deterministic midpoint diagnostic.  Current GRHayL NRPyLeakage consumers
interpolate only `muhat,mu_e,mu_p,mu_n,Xn,Xp`; qualification tests those six
values and downstream range errors directly.

## Verification

Routine tests use an asymmetric official-schema-shaped synthetic fixture and
require 100% production Python line and branch coverage:

```sh
compose-venv/bin/python -m pip install coverage
compose-venv/bin/python -m coverage run --branch --source=tools/compose \
  --omit='Unit_Tests/*' -m unittest discover \
  -s Unit_Tests/compose -p 'test_*.py'
compose-venv/bin/python -m coverage report --fail-under=100
compose-venv/bin/python -m coverage xml -o compose-coverage.xml
```

After producing a table, build GRHayL with HDF5 and run the auto-discovered C
integration test:

```sh
./configure -r
make test/unit_test_tabulated_eos_compose
LD_LIBRARY_PATH="$PWD/build/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" \
  ./test/unit_test_tabulated_eos_compose TABLE.h5
```

The C test independently reads every output node and all 19 fields, verifies
the unchanged loader mapping and units, derived enthalpy, relativistic
sound-speed path, off-grid storage-space interpolation, `eps/P/S/h`
temperature inversions from distinct valid initial guesses at endpoints and
an interior point, the six-value
NRPyLeakage callback order and range failures, and memory cleanup.  It also
checks a no-fallback Con2Prim recovery, analytic tabulated flux and source
goldens, and all eight combined-NRPyLeakage outputs against fixed regression
goldens for the two qualified table dimensions.

On the qualified table-141 artifact, regularization changed at most
0.0341587 in `logenergy`, 0.0338011 in `logpress`, and 0.00547615 in a mass
fraction; maximum entropy change was `1.1608683721e10`.  The deterministic
energy shift was `30.0518 MeV` per baryon.  Regularization applied the high
causal bulk bound at 221133 nodes and the low bound at 1578 nodes, created
2184728 absent-heavy sentinels, and clipped 257 heavy-charge roundoff values
by at most `6.16751e-12`.  Maximum midpoint off-grid charge residual was
`0.0408931`.  These material changes are why the output has a distinct profile
name and must not be represented as the official raw EOS.
