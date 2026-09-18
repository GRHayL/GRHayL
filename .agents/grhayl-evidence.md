# GRHayL Evidence And Validation

Use applicable sections with [the protocol](review-protocol.md), not a blanket
checklist. Follow [KB Index](../wiki/index.md), [Catalog](../wiki/catalog.md), and
exact owners; this reference adds no workflow or budget.

## Authority And Interfaces

Source, headers, recursive `make.code.defn`, `configure`, dispatch, tests, and
invoked binaries establish behavior; [derivation](../docs/raw/derivation.md) and
owning Doxygen sources supply equation/variable intent. KB pages route/synthesize;
[Contradictions](../wiki/contradictions.md) identifies unresolved seams. Votes cannot
choose maintainer intent.

Keep declaration, definition, manifest/configure selection, installation, dispatch,
test selection, compile/link, execution, CI selection, and downstream proof distinct.
Use [Public API Map](../wiki/public-api-map.md) and [Source Map](../wiki/source-map.md)
for affected guards, consumers, manifests/install lists, docs/tests, and downstream
contracts.

## Builds And Numerical Tests

`make tests`/`make datagen` compile/link; invoke binaries separately. Distinguish
generator execution, fixture production, replay, perturbation, checked returns,
and assertion strength. Process success or self-consistency is not numerical proof.
Justify oracles/tolerances; do not loosen tolerances or regenerate references just
to pass. For HDF5/source-selection changes, inspect `configure`, `scripts/parser`,
recursive manifests, and generated `SRC`, `TEXES`, `DGEXES`, `IHDS`; check applicable
default and `--disable-hdf5` modes. Not all tabulated code necessarily disappears.

Run [.github/run_tests.sh](../.github/run_tests.sh) only in a disposable checkout
after inspecting downloads, cleanup globs, and failure leftovers. While
[Contradictions](../wiki/contradictions.md) marks `generate_makefile.sh` output broken,
inspect it only in isolation; never run `make` from malformed output.

## Generated Code And Doxygen

Follow [Generated Boundaries](../wiki/generated-boundaries.md). Use a unique temporary
Doxygen `OUTPUT_DIRECTORY` with correct repo-relative inputs; never generate beside
`docs/raw/`. Compare exit status/warnings using matching versions/configurations.
For Flux_Source/NRPyLeakage derived C, establish supported generator routes first.
Checked-in C/headers own behavior; never invent commands or claim reproducibility
across documented drift/unknown entry points. Include authorized generator/product
companions; avoid unrelated regeneration.

## CI And Downstream

Workflow/action presence proves selection, not historical execution. Consider
relevant events, path filters, compiler, OS, and HDF5 variants. Upstream/ET_Legacy
tests or static GRHayLib parity do not prove Cactus build, schedule, parameters,
lifecycle, or runtime. Follow
[downstream coordination](../wiki/workflows.md#downstream-grhaylib-impact) before
editing `implementations/GRHayLib/`; report unavailable environments using the
[verification boundary](../wiki/implementations/grhaylib/verification-and-drift.md).

## KB And Agent Instructions

Prepare one candidate, not competing KBs. Follow [KB contracts](../wiki/index.md#page-contract)
and [KB Checks](../wiki/lint/CHECKS.md), which are manual examples, not a canonical
linter. Reopen exact ground truth for changed claims; update affected routes, links,
owner pages, catalog/source/test/API maps, and contradictions only as needed.
KB-only scope does not authorize source, Doxygen, workflow, or downstream repairs;
retain path-specific read-only boundaries, including Reconstruction documentation.

Check agent metadata, links/anchors, activation, precedence, coverage, budgets, and
exit/delivery paths. Use affected link/policy checks and `git diff --check` when Git
is available, not unrelated numerical suites. Keep repo-relative KB links and the
runtime-path exception in [AGENTS.md](../AGENTS.md#agent-execution-contract).
No source-tracking checksums, hashes, digests, VCS revision pins, file or source
counts, `mtime`, stored fingerprints, date stamps or timestamps as KB metadata,
or separate KB logs. Technical hash facts may be documented only
as reviewed domain facts, never stored digest values. Operational checkpoints
are task state, not permanent KB logs. Do not import NRPy-only tools, schema, or
file protections.
