# NRPyLeakage follow-up fix plan

## Objective

Close the actionable remainder established by the two independent audits of
`issuesCC.md`, without reopening claims that the evidence rejected or
overstated. The implementation should be minimal, mutation-sensitive, and keep
the existing NRPyLeakage public ABI.

The required work is:

1. define an honest single-species endpoint contract (C-3);
2. reject floating-point build modes that defeat nonfinite handling (C-6);
3. add the missing independent public diffusion oracle (C-1/C-4);
4. make the classifier, sphere-caller, generator-pairing, and EOS-failure
   regressions persistent (C-5, C-8, C-9, C-10);
5. complete the installed API, runner, test, and licensing documentation
   (C-13, C-15, and the valid C-16 runner subclaim).

Minimality rule: prefer a local condition, an existing error code, and an
existing test executable. Add no public ABI, new production abstraction, new
test executable, fixture-format change, or broad refactor unless a focused
negative control proves it is necessary. Tests are authorized to land in the
main branch, but should be added to the existing NRPyLeakage test sources.

## Fixed scope and non-goals

### In scope

| Finding | Accepted problem | Planned resolution |
| --- | --- | --- |
| C-3 | A normalized single-species endpoint returns success while both charged-current orientations are zeroed. | Define and test the exact endpoint as the analytic continuous extension of the installed density-derived closure: zero charged-current products, occupied-species scattering retained, and no divergent reaction shift formed. |
| C-6 | Accepted compiler flags can remove the exceptional-value semantics required by the nonfinite-output contract. | Reject known unsafe flags and detect the effective compiler mode in the supported `configure` route. |
| C-1/C-4 | Helper identities do not independently validate public heavy-species diffusion wiring, despite documentation claiming they do. | Reconstruct the public thick-limit result from public opacity, equilibrium energy density, and optical depth. |
| C-5 | Runtime values do not prevent reintroduction of effective-type aliasing, and the portable fallback matrix is incomplete. | Add a strict-aliasing compile regression and expand portable classifier cases. |
| C-8 | The corrected sphere caller remains insensitive to neighbor reversal. | Exercise the sphere's actual gather/call setup with asymmetric data and document stencil layout. |
| C-9 | Paired inputs are neither asserted nor protected at EOS boundaries. | Preserve base/perturbed vectors, bound endpoint perturbations, and assert the pairing invariant. |
| C-10 | The five repaired EOS checks lack failure injection; the initial lookup can truncate trusted output before it is validated. | Propagate test-local errors, inject all five failures, and perform the initial checked lookup before opening the fixture. |
| C-13 | Installed declarations omit error and output-validity contracts. | Add complete Doxygen blocks to the three public EOS-dependent routines. |
| C-15 | Licensing summaries omit or misdescribe the FDINT BSD-3 exception. | Correct the guide, catalog route, and GRHayLib README. |
| C-16, runner subclaim | The documented local command transcript omits the fallback executable and blurs local versus workflow execution. | Correct the tests-and-fixtures page. |

### Explicitly out of scope

- Do not clamp, floor, or otherwise change the reaction shift `q` for C-2.
  The detailed-balance defect claimed there was not established. Existing
  dense-matter accuracy limitations remain a scientific qualification.
- Do not change the free-nucleon scattering model (C-7), fixture tolerances
  (C-11), or the legacy Fermi check (C-12).
- Do not add an EOS mean-field/effective-mass interface. A physically usable
  single-species limiting model needs such a separately designed extension.
- Do not add `warn_unused_result`, merge comparison wrappers, rename private
  helpers, hide private Doxygen members, or add an unsupported Fermi-family
  cross-check.
- Do not make no-HDF5 or SRO141 a new CI matrix requirement in this change.
  They remain useful validation/manual routes, not demonstrated current CI
  failures.
- Do not regenerate or publish binary fixtures unless the focused work proves
  that their intended values must change. Publication is separate authority.
- Do not add a new test executable or fixture format. Extend the existing
  NRPyLeakage physics, fallback, sphere, luminosity, and optically thin tests.
- Do not change Cactus scheduling, CCL interfaces, or the GRHayLib source
  wiring.

## Design decisions

### 1. Single-species endpoints use the closure's analytic product limit

The current EOS callback supplies free fractions and chemical potentials, but
not enough mean-field information to define a separate interacting-matter
endpoint model. It does, however, define a unique continuous endpoint for the
installed density-derived closure. As the minority fraction tends to zero,
the kinetic degeneracy difference diverges logarithmically, the reverse
transition overlap vanishes exponentially, and every overlap-times-shifted-
beta-moment product tends to zero. The occupied-species scattering population
and the non-beta channels remain finite.

Use this contract after the helper's existing roundoff normalization:

- `X_n > 0 && X_p > 0`: retain current behavior exactly;
- `X_n == 0 && X_p == 0`: return success with no free-nucleon charged-current
  targets; pair/plasmon processes remain available;
- exactly one of `X_n`, `X_p` is zero: return success with the exact overlap
  limit, retain occupied-species scattering, and apply the analytic zero limit
  of the charged-current products without constructing `q`;
- a small accepted negative fraction that normalizes to zero follows the same
  rule;
- every strictly positive fraction remains accepted. Add no numerical floor.

Centralize the overlap/scattering rule in
`GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_nucleon_blocking.h`, inside
`NRPyLeakage_compute_nucleon_blocking()`. Keep the three callers' both-positive
beta-moment guards and document them as the analytic product-limit branch. Add
no new enum, signature, or repeated numerical policy.

Do not describe this closure limit as proof of general interacting-EOS endpoint
accuracy. Do not invent a finite `q`, clamp fractions, or floor compositions.

### 2. Unsupported floating-point semantics fail at build time

NRPyLeakage's nonfinite status cannot be recovered after compiler
transformations have converted exceptional input into apparently finite
arithmetic. Configuration must therefore reject the build instead of trying to
repair its outputs.

Use two complementary checks:

1. In `configure`, reject known unsafe user tokens: `-ffast-math`, `-Ofast`,
   `-ffinite-math-only`, `-funsafe-math-optimizations`, `-fassociative-math`,
   `-freciprocal-math`, `-fno-signed-zeros`, and `-fno-trapping-math`. Reject
   umbrella modes even when later flags appear to negate only one of their
   consequences; the project does not claim that a partial reversal restores
   all required semantics. Keep this list centralized so configure help,
   diagnostics, and tests cannot drift.
2. After compiler selection and final flag assembly, compile a small probe that
   errors when `__FAST_MATH__` is defined or `__FINITE_MATH_ONLY__ > 0`. This
   catches effective modes and compiler aliases rather than relying solely on
   spelling.

Document the restriction in `configure -h`, `wiki/build-and-ci.md`, and the
NRPyLeakage implementation/API documentation. Do not state that appending only
`-fno-finite-math-only` makes arbitrary fast-math combinations supported.
Direct Cactus builds do not use `configure`; document that they must apply the
same flag restriction, but do not add a second enforcement mechanism without a
reproduced downstream failure.

## Implementation work packages

### WP1 — Endpoint contract and public regression matrix

Files:

- `GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_nucleon_blocking.h`
- `Unit_Tests/unit_test_nrpyleakage_physics.c`

Steps:

1. Preserve the exact overlap limits after existing fraction validation and
   roundoff normalization: occupied-to-empty equals the occupied fraction and
   the reverse overlap is zero.
2. Preserve the both-zero successful return, occupied-species scattering, and
   the existing strictly-positive calculation.
3. Parameterize the test EOS callback so tests can supply independent `X_n`
   and `X_p` values.
4. Exercise both single-species orientations through all three public APIs:
   opacities, luminosities, and the combined opacity/source routine.
5. Cover exact zero and a permitted roundoff-negative value normalized to zero.
6. Require success and finite outputs; explicitly require every standalone and
   combined opacity moment to retain scattering above the fallback floor and
   require non-beta emission to remain active.
7. Add a both-zero public case requiring success and finite outputs with zero
   charged-current contribution.
8. Trace representative tiny-positive cases down to `1e-300` in both
   orientations and require convergence to the exact public endpoint through
   all three APIs. This proves a continuous extension rather than an
   undocumented threshold.
9. Update the old direct-helper endpoint assertions to the same contract.

Before accepting this behavior, run all available SLy4 replays and inspect
whether supported fixture states reach a normalized single-species endpoint.
They do; therefore reject the incompatible hard-error proposal and require all
unchanged replays to pass with the analytic endpoint extension. Do not hide the
endpoint with a fraction floor or rebaseline fixtures around it.

Acceptance:

- exactly-one-zero and normalized-to-zero inputs agree with the analytic
  tiny-positive limit through all public paths;
- occupied-species scattering, non-beta processes, both-zero, and strictly-
  positive behavior remain supported;
- no existing positive-composition golden changes.

### WP2 — Independent public heavy-species diffusion oracle

File:

- `Unit_Tests/unit_test_nrpyleakage_physics.c`

Build an oracle from observable public quantities rather than internal emission
helpers or regenerated constants:

1. Use a mock state whose electron and heavy-species free rates differ.
2. Call the public opacity routine and obtain the heavy-species energy opacity
   `kappa.nux[1]` in geometric units.
3. Keep number optical depths at zero. Obtain the free public heavy luminosity
   at zero energy depth, divide out
   `NRPyLeakage_units_cgs_to_geom_Q * W * alpha^2 * sqrt(det(gamma))` to recover
   the single-species cgs `Q_free`, then evaluate a large heavy energy depth
   such as `tau_E = 1000`.
4. Independently compute

   \[
   U_{\mathrm{nux}} =
   \frac{4\pi T^4 F_3(0)}{(hc)^3},
   \qquad
   Q_{\mathrm{diff}} =
   \frac{\kappa_{\mathrm{geom}}}{L_{\mathrm{unit}}}
   U_{\mathrm{nux}}\frac{c}{6\tau_E^2}.
   \]

   Obtain `F_3(0)` through the public Fermi-integral routine. Keep every unit
   conversion explicit in the test.
5. Reconstruct the finite-depth effective rate, not merely its asymptote:

   \[
   Q_{\mathrm{eff}} =
   \frac{Q_{\mathrm{free}}}
        {1 + Q_{\mathrm{free}}/Q_{\mathrm{diff}}}.
   \]

6. Convert this result back with the same public
   `NRPyLeakage_units_cgs_to_geom_Q * W * alpha^2 * sqrt(det(gamma))` factor and
   compare it with `lum.nux`.
7. Repeat through the combined public routine. Isolate the change caused by
   changing only the heavy energy depth and include the factor of four heavy
   species in the matter energy source.
8. Derive a tight rounding tolerance from the explicit reconstruction. Do not
   use a loose asymptotic tolerance or add another implementation-derived
   golden.

Required mutation proof during implementation: restore the wrong
electron-species denominator at each heavy call site. The new public oracle
must fail even when helper identities remain unchanged. Record the mutation
result in the change review; do not retain mutation code.

Acceptance:

- both public emission paths agree with the independent finite-depth
  reconstruction;
- swapping the heavy denominator species fails by a clear margin;
- the test separately protects species selection, opacity conversion, energy
  depth, and the heavy-species multiplicity.

### WP3 — Build-semantics enforcement

Files:

- `configure`
- `.github/workflows/github-actions-Ubuntu-gcc.yml`
- `wiki/build-and-ci.md`
- `wiki/gems/neutrinos/api-and-data.md`
- `wiki/gems/neutrinos/implementation-flow.md`

Steps:

1. Add the token and effective-mode checks described above.
2. Ensure the diagnostic names NRPyLeakage's exceptional-value/nonfinite-output
   requirement and the offending effective mode.
3. Add one focused Ubuntu GCC workflow step, on one matrix leg only, that
   requires configuration failure for each of:
   `--cflags=-ffast-math`, `--cflags=-Ofast`, and
   `--cflags=-ffinite-math-only`, plus the representative individual transform
   `--cflags=-funsafe-math-optimizations`.
4. Require default configuration to continue succeeding. Also test an ordinary
   custom flag so the check does not reject all user flags.
5. Run every negative configuration in its own clean temporary checkout/build
   and match the intended unsupported-floating-mode diagnostic. A failure due
   to an existing `build/.check`, missing dependency, or unrelated configure
   error must fail the workflow check rather than count as success. Run the
   positive control in another clean checkout.

Do not solve this by merely moving `-fno-finite-math-only` after user flags.

Acceptance:

- the listed unsupported modes and macro-detectable aliases cannot produce an
  installed library through `configure`;
- direct-build documentation states the same unsupported flag boundary without
  claiming that `configure` enforces Cactus compilation.

### WP4 — Classifier regression coverage

Files:

- `Unit_Tests/unit_test_nrpyleakage_classifier_fallback.c`
- `Unit_Tests/unit_test_nrpyleakage_physics.c`
- `GRHayL/include/ghl_nrpyleakage.h`
- `.github/workflows/github-actions-Ubuntu-gcc.yml`

Steps:

1. Expand the portable fallback executable to test positive and negative zero,
   ordinary finite values, maximum finite values, positive and negative
   subnormals, both infinities, and NaN signs.
2. Do not promise signaling-NaN representation behavior on the portable
   fallback path.
3. Add an explicitly test-only
   `GHL_NRPYLEAKAGE_FORCE_PORTABLE_CLASSIFIERS` switch around the binary64
   classifier selection in `ghl_nrpyleakage.h`. Define it before including the
   header in the fallback executable. This avoids relying on whether a
   conforming implementation redefines `UINT64_C` during nested header
   inclusion. Do not present the old `#undef UINT64_C` as undefined behavior;
   it is being replaced for deterministic branch selection.
4. Reuse `unit_test_nrpyleakage_physics.c` as the strict-aliasing compile
   translation unit; it already instantiates both classifiers. Compile it with
   the configured include paths on one Ubuntu GCC route, without linking a
   second executable, using at least:

   ```text
   -O2 -fstrict-aliasing -Wstrict-aliasing=2 -Werror=strict-aliasing
   ```

5. Confirm as a temporary negative control that restoring the old pointer puns
   fails this compile while the current `memcpy` implementation passes.

Acceptance:

- runtime matrices cover both the binary64 and portable classifier paths;
- the persistent compile check rejects the original effective-type violation;
- no project-wide warning policy is changed.

### WP5 — Sphere caller wiring regression

Files:

- `Unit_Tests/unit_test_nrpyleakage_constant_density_sphere.c`
- `GRHayL/include/ghl_nrpyleakage.h`

Refactor only the test-local caller setup:

1. Extract a static helper that computes the center/six neighbor indices,
   gathers opacity and depth values from the sphere's arrays, and makes the
   production-like `NRPyLeakage_optical_depths_PathOfLeastResistance()` call.
2. Use that same helper in the existing physical sphere loop.
3. Add a small focused case with distinct opposite-face opacities/depths and
   unequal `[minus, center, plus]` metric values in every direction.
4. Hand-compute the six candidate path costs and assert the expected minimum.
5. Temporarily reverse one gathered neighbor pair and require the focused case
   to fail. Do not alter the physical flat-sphere fixture merely to create this
   sensitivity.

Add public Doxygen stating:

- `dxx` is `[dx, dy, dz]`;
- each metric stencil is `[minus face, center, plus face]`;
- each minus/plus opacity must stay paired with the depth from the same face.

Acceptance: the actual sphere gather/call path, not only a direct kernel call,
fails a re-reversal mutation.

### WP6 — Bounded, explicitly paired luminosity generation

File:

- `Unit_Tests/unit_test_nrpyleakage_luminosities.c`

Steps:

1. Preserve the existing base values in `const` locals and compute separate
   perturbed locals. Do not introduce a new input structure or mutate the base
   values in place.
2. For `rho`, `Y_e`, and `T`, use a bounded relative-perturbation helper:
   try the sampled sign; if it exits the closed EOS domain, reflect the sign;
   fail loudly if neither candidate is in range.
3. For all components, assert finiteness and the componentwise invariant

   ```text
   abs(perturbed/base - 1) <= 1e-14 + rounding_allowance
   ```

   with a separate exact-zero rule that preserves zero. Base values in this
   generator are normally nonzero, but the helper contract must still be
   explicit.
4. Assert `rho`, `Y_e`, and `T` remain inside the EOS bounds immediately before
   the public luminosity call.
5. Add direct generic-helper cases at its lower and upper bounds with both
   perturbation signs; do not duplicate the same helper test for each field.
6. Preserve the fixture layout. The input invariant is checked during
   generation; it need not be encoded into the binary replay format.

Sign reflection is selected over narrowing the sampling interval or clamping:
it handles legal RNG endpoints without creating clamped duplicate values and
preserves current fixed-seed bytes when no draw leaves the domain. It does force
an inward perturbation at an exact bound; document that deliberate boundary
bias rather than claiming none exists.

Acceptance:

- every written output pair has a checked matching input pair;
- legal RNG endpoints cannot send EOS-bounded inputs out of range;
- scratch regeneration is byte-identical to existing fixtures unless an
  independently reviewed intentional change is found.

### WP7 — EOS failure injection and non-destructive initial validation

File:

- `Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c`

Steps:

1. Change test-local `rhs()`, `rk4_step_ode()`, and one-fixture generation
   helpers to return `ghl_error_codes_t`. Keep `ghl_abort_if_error()` only at
   the outer executable boundary for ordinary generation/replay.
2. Introduce test-local callback indirection/counters for the initial
   `compute_eps_from_T` lookup and the four RK `compute_T_from_eps` lookups.
3. Inject failure at RK stages 1 through 4. Require the original error and
   exactly `0`, `1`, `2`, and `3` completed RHS calls, respectively. The failed
   stage must never evaluate its RHS.
4. Inject failure at the initial energy lookup. Require unchanged output
   sentinels and no dependent record write.
5. Move the initial checked `compute_eps_from_T` lookup before `fopen(...,
   "wb")`. Open the output only after that lookup succeeds. This is the smallest
   fix for the reproduced zero-byte truncation; do not add temporary-file,
   backup, or two-file transaction machinery in this change.
6. Add a scratch-directory regression containing a marker canonical file.
   Inject the initial lookup failure and require the marker to remain
   byte-identical.

Later RK failures may still leave a partially generated key-0 artifact after
the file is intentionally opened. That is not trusted published data and was
not the reproduced failure-atomicity defect. A fully transactional generator is
a separate hardening change, not required here.

Acceptance:

- all five repaired lookup sites have site-specific persistent negative tests;
- no value derived from a failed lookup is written;
- an initial lookup failure preserves an existing canonical fixture.

### WP8 — Installed API and model documentation

Files:

- `GRHayL/include/ghl_nrpyleakage.h`
- `wiki/gems/neutrinos/api-and-data.md`
- `wiki/gems/neutrinos/implementation-flow.md`
- `wiki/gems/neutrinos/physics-and-eos-contract.md`
- `wiki/gems/neutrinos/nucleon-blocking-and-eos-conventions.md`

Add full Doxygen blocks to the three public EOS-dependent routines. For each,
document:

- parameter units and output meaning;
- success and propagated EOS/Fermi/HDF5 errors;
- `ghl_error_nrpyleakage_blocking`, plus the successful exactly-one-zero
  analytic-limit rule;
- which early failures occur before any output write;
- `ghl_error_nrpyleakage_nonfinite_output` and its deterministic finite
  fallback values;
- that sanitized fallback values are not a certification of physical validity;
- opacity floors versus zero source/luminosity fallbacks.

Keep fatality/recovery policy caller-owned. Do not describe the nonfinite status
as merely advisory.

Update the KB pages to explain why exact single-species populations retain the
analytic zero limit of the closure's charged-current products while scattering
and non-beta channels remain active. Preserve the documented common-bare-mass
approximation and dense-matter limitations; make no general physical-accuracy
claim.

### WP9 — Test/runner documentation

Files:

- `wiki/gems/neutrinos/tests-and-fixtures.md`
- `wiki/test-map.md`
- `wiki/tests/unit-test-coverage-and-gap-matrix.md` if its coverage summary is
  changed by the implementation

Changes:

1. Attribute caller species/unit protection to the new public thick-limit
   reconstruction. Describe helper identities only as formula-level checks.
2. Describe `.github/run_tests.sh` as a broad local replay driver. Workflows
   invoke the relevant neutrino executables directly; they do not call that
   script.
3. Add `unit_test_nrpyleakage_classifier_fallback` to the documented local
   command transcript.
4. Add the compile-only classifier check, sphere caller regression, endpoint
   matrix, pairing invariant, and EOS-failure injections to the test map.
5. Distinguish a configured command from an observed CI run.

Do not edit `.github/run_tests.sh` for the fallback test: it already invokes
that executable.

### WP10 — Licensing summaries

Files:

- `docs/raw/license.md`
- `wiki/catalog.md`
- `implementations/GRHayLib/README`

Changes:

1. Add the FDINT-derived minimax routines as BSD-3-Clause material in the main
   license guide and route readers to the exact source notice and root
   `THIRD_PARTY_NOTICES`.
2. Add `THIRD_PARTY_NOTICES` to the catalog's licensing ground-truth route.
3. Replace the GRHayLib README's blanket “BSD2 for all other code” wording with
   the actual compiled-code exceptions: GPL for Noble-derived code, CC BY-SA
   for Palenzuela/Newman-derived code, Boost Software License 1.0 for `toms748`,
   BSD-3-Clause for WENO-Z and FDINT-derived code, and BSD-2-Clause for the
   remaining GRHayL code.

Do not duplicate complete license texts. The source-local headers and installed
root notice remain authoritative. No source or install-rule change is expected.

## Implementation order

1. Add focused failing tests for C-3, C-1/C-4, C-5, and C-8.
2. Make the endpoint-limit contract explicit and implement the build-semantics guards.
3. Prove the public-oracle, pointer-pun, and neighbor-reversal negative
   controls; remove all temporary mutations.
4. Make the focused luminosity and optically thin test-generator changes.
5. Update installed API comments, KB pages, runner text, and licensing text.
6. Run the complete validation matrix below in disposable checkouts/builds.
7. Reinspect the final diff for accidental fixture, generated-file, or
   unrelated changes.

This ordering makes semantic failures visible before documentation is rewritten
and separates endpoint-contract verification from test-generator hardening.

## Validation matrix

Use disposable checkouts for configuration/build variants and any command that
generates or downloads fixtures. The repository runner removes broad root-level
file globs and must not be run over user-owned artifacts.

### Static and configuration checks

- `git diff --check`
- `sh -n configure .github/run_tests.sh`
- default GCC and Clang configuration succeeds;
- ordinary custom warning/optimization flags succeed;
- every listed unsupported flag fails in an isolated configure tree with the
  intended diagnostic, including the representative individual-transform CI
  case `--cflags=-funsafe-math-optimizations`;
- parse changed workflow YAML and run `actionlint` when available.

### Build and runtime checks

Default HDF5 configuration:

1. `./configure -r`
2. `make -j<N> grhayl tests datagen`
3. run `unit_test_nrpyleakage_physics`;
4. run `unit_test_nrpyleakage_classifier_fallback`;
5. run the luminosity, optically thin gas, and constant-density sphere SLy4
   replays when their published fixtures/table are available;
6. run the relevant expected-error cases for disabled HDF5 and invalid blocking inputs.

No-HDF5 configuration:

1. configure with `--disable-hdf5`;
2. build table-free tests;
3. run the physics and classifier-fallback executables.

### Focused regression checks

- all exact/roundoff endpoint cases across all three public APIs;
- both-zero and tiny-positive boundary cases;
- independent public heavy finite-depth reconstruction in luminosity and
  combined-source paths;
- wrong-species denominator mutation fails the new oracle;
- strict-alias compile passes current code and fails restored pointer puns;
- portable classifier matrix covers signed zero, subnormals, infinities, and
  NaNs;
- sphere asymmetric setup passes, and a reversed neighbor mutation fails;
- forced legal RNG endpoints remain paired and EOS-bounded;
- all five EOS failure injections return before the prohibited RHS/write;
- a marker fixture survives injected initial-lookup failure unchanged.

### Fixture checks

- run both key-0 generators in a disposable directory;
- generate twice from the same seed and compare the two runs;
- compare scratch results with available published fixtures;
- expect no fixture update from this plan;
- if bytes differ, identify the first differing input/output and review the
  cause before seeking separate publication authority. Never silently
  rebaseline.

### Documentation, install, and downstream checks

- run Doxygen with a unique temporary `OUTPUT_DIRECTORY`; inspect warnings and
  confirm the three public declarations and stencil contract render;
- install into a temporary prefix and verify the public header plus
  `THIRD_PARTY_NOTICES`;
- compile a small consumer against the installed header;
- verify GRHayLib's source symlinks still resolve and its README routes to the
  central notices;
- no Cactus runtime claim is made without an actual Cactus build/run.

## Delivery gates

The change is ready only when all of these are true:

- the endpoint contract is explicit, uniformly returned, and tested at public
  boundaries;
- supported table replays do not expose an unresolved endpoint/model conflict;
- all compiler modes declared unsupported by this plan fail before producing a
  library through `configure`, and direct-build limitations are documented;
- each repaired historical defect has a mutation-sensitive or failure-injected
  persistent check at the affected integration layer;
- the initial lookup cannot truncate a trusted fixture, and generated
  luminosity inputs are bounded matched pairs;
- installed API docs accurately distinguish untouched outputs, sanitized
  fallbacks, and physical success;
- licensing and runner summaries match repository ground truth;
- no rejected issue, unrelated refactor, generated Doxygen output, or fixture
  rebaseline enters the diff.

## Known limits after completion

- The change will not establish interacting-EOS accuracy at large finite
  `|q/T|`; it will preserve and accurately document that qualification.
- Exact single-species states use the continuous endpoint of the installed
  closure; this does not establish the endpoint rates of a richer interacting
  EOS model.
- Full GitHub Actions observation, external DD2/BNS_NURATES qualification,
  SRO141 manual comparison, non-LP64 execution, and Cactus runtime testing
  remain unperformed unless separately scheduled.
- Later key-0 RK or process failures can still leave a partial generation
  artifact after opening; full transactional publication is intentionally
  outside this minimal fix.
