---
name: trialectic
description: Use when the user explicitly asks to use tri, engage tri, run tri, run or engage the GRHayL trialectic, requests three-seat independent GRHayL review, or repository policy requires it. Produce or review one GRHayL change through three fresh-context seats, repo-local evidence, proportionate validation, at most three candidate cycles, and unanimous acceptance except for one tightly bounded terminal correction of mechanical defects. Do not invoke for routine edits without an explicit trigger or policy requirement.
---

# GRHayL Trialectic

Use three independent seats to produce or review one GRHayL candidate. Keep initial work independent, then resolve findings against the frozen brief, claim-specific repository ground truth, and executable checks. Agreement is never stronger than evidence.

Root is lead. Root owns discovery, scope, synthesis, the canonical candidate, validation, installation, and reporting. Delegated seats write only in separate temporary areas and never to live targets.

## Invariants

- Preserve the user's request, authorization boundaries, unrelated work, and concurrent changes.
- Use exactly three fresh-context seats in every delegated phase. If three complete independent results are unavailable, stop rather than degrade silently.
- Root launches every wave. Neither root nor seats recursively invoke Dialectic or Trialectic; seats do not create subagents.
- Give all seats the same frozen brief, evidence boundaries, acceptance criteria, candidate, and validation record. Vary only their complementary roles. Do not reveal peer work before all initial results freeze.
- Root resolves findings against evidence. Do not add peer debate, cross-critique, rebuttal, or vote-seeking waves.
- Permit `NO ACCEPTABLE CANDIDATE`. Never weaken requirements, force synthesis, treat silence as approval, or resample a substantive block as agent failure.
- Treat applicable `AGENTS.md` files as ambient instructions. Modify one only when the user explicitly targets it and governing instructions permit the change.
- Never overwrite concurrent work. Material post-freeze drift ends the invocation.

## GRHayL evidence model

Start at `wiki/index.md`, then use `wiki/catalog.md`. Read the relevant owner leaf and direct ground-truth files. Use `wiki/source-map.md`, `wiki/workflows.md`, `wiki/test-map.md`, `wiki/public-api-map.md`, `wiki/generated-boundaries.md`, `wiki/change-impact.md`, and `wiki/contradictions.md` when their boundaries apply.

Authority is claim-specific. Current executable behavior comes from the exact source, public headers, recursive `make.code.defn` manifests, `configure`, dispatch or initialization, tests, and invoked binaries. Equations and variable definitions may require `docs/raw/derivation.md` or owning Doxygen source. CI selection comes from workflows and composite actions. Downstream claims require `implementations/GRHayLib/` and, for build or runtime proof, an appropriate Cactus environment. KB pages route and synthesize; they do not resolve conflicts alone.

Keep support layers distinct: declaration, definition, manifest inclusion, configure selection, installed exposure, dispatch or initialization, test selection, compilation or link, exact execution, CI selection, and downstream runtime proof. Presence in one layer proves no other layer.

Consult `wiki/contradictions.md` before resolving an unsafe seam. State observable behavior and limitations; agreement cannot choose unresolved maintainer intent.

## Choose mode

Record the least expensive useful mode before delegation.

- `review` — Default. Root prepares one complete candidate; three seats review it. Use when repository patterns and validators substantially determine the change.
- `design` — Three seats independently propose complete designs; root selects or synthesizes one coherent design, prepares a candidate, then uses a fresh triad for review. Use only when materially different valid designs remain after discovery.
- `tri-build` — Three seats independently build complete isolated candidates; root selects or synthesizes one coherent candidate, then uses a fresh triad for review. Use only when independent implementations provide a meaningful oracle or reduce material scientific or integration risk.

Documentation, KB, configuration, routine fixes, and established source patterns normally use `review`. An explicitly requested competing information architecture may use `design`, but seats never choose source authority by vote.

## Budget and stop rules

One invocation has one frozen scope round. After freeze, a material change to targets, intent, acceptance criteria, authoritative inputs, repository rules, generated ownership, or validation obligations requires `RESTART REQUIRED`.

Limits: one proposal or tri-build wave; three candidate-review cycles; one review-and-decision wave per cycle; one root correction batch between consecutive cycles; one agent-recovery wave for mechanical or malformed output; and one retry for a diagnosed transient validation failure.

Two cycles are the normal ceiling. Open cycle three only after concrete progress and with one bounded correction or authorized evidence acquisition plus a deciding check. Stop when a blocker or evidence gap recurs without materially new evidence, a deterministic failure has no bounded fix, required evidence is unavailable, scope must expand, candidate identity cannot bind, drift occurs, or the next wave exceeds budget.

Run each planned check once per candidate identity. Bounded read-only diagnosis may refine a failing command, but never rerun an unchanged deterministic failure to seek a pass. On correction, rerun every check whose inputs may intersect changed bytes.

## Freeze brief and candidate

Frozen brief states:

- request, acceptance criteria, explicit non-goals, and exact create/change/delete targets;
- relevant live baselines, preservation rules, assumptions, and compatibility requirements;
- bounded source, header, build, test, documentation, CI, generated, KB, and downstream evidence;
- applicable contradictions and generated-output ownership;
- selected mode and complementary roles;
- validation commands, expected evidence, unavailable checks, untested environments, and tolerance rationale; and
- one root candidate area plus one separate temporary area per seat, outside live targets.

A seat needing material evidence outside the brief returns `SCOPE EXPANSION REQUEST` with exact path or bounded area, reason, and impact. If necessary, root stops with `SCOPE INCOMPLETE`; root does not expand frozen scope mid-run.

Before review, root completes formatting, fixers, and authorized regeneration; freezes one canonical snapshot; inventories all operations; binds validation and decisions to its identity; and prevents validators from mutating it. Run mutating tools on expendable identical copies. Any candidate mutation or substantive change to decision evidence invalidates prior decisions and normally consumes the next cycle with fresh reviewers.

For ordinary files, identity may use immutable copies, direct byte equality, or content digests. For `AGENTS.md`, `wiki/**`, and future GRHayL KB/governance sources, follow repository policy: no source-tracking hashes or hashing of sources, no `mtime`, no stored fingerprints, no maintenance log, and no unnecessary dates. Use frozen copies or direct comparison; necessary dates use `MM-DD-YYYY`.

## Complementary seats

### Seat 1: physics, numerics, and behavior

Check applicable GRMHD, EOS, neutrino, reconstruction, and induction assumptions; equations and variable conventions; units, indices, signs, densitization, face or stagger ownership, bounds, failure paths, numerical stability, tolerances, and scientific oracles. Trace behavior through exact sources and tests. Require justified equivalence, reference comparison, convergence, or focused runtime checks only when relevant.

### Seat 2: API, build, and integration

Check public headers, structs and enums, declarations versus definitions, manifests, `configure` selection, HDF5 guards, dispatch, installation surface, API/ABI compatibility, ownership boundaries, direct consumers, and impact across Core and gems. Prefer established GRHayL patterns and the smallest sufficient implementation.

### Seat 3: tests, docs, CI, generated, KB, and downstream audit

Trace every acceptance criterion to repository evidence. Check targeted and error tests, fixture lifecycle and assertion strength, Doxygen and KB routes, CI selection, generated provenance and reproducibility limits, stale references, packaging, and `implementations/GRHayLib/` CCL/build/lifecycle impact. Reject support or release claims that outrun evidence and expose omissions missed by local correctness review.

For non-code work, use a first-principles analyst, a domain-and-implementation expert, and an adversarial repository-evidence editor.

## GRHayL validation rules

Choose checks proportional to changed paths and risk. Record working directory, command, exit status, decisive output, candidate identity, unavailable dependencies, untested modes, and tolerance rationale.

- Isolate commands whose output or cleanup could contaminate live targets or the canonical snapshot. Use disposable trees for full builds/runs, fixture generation, installs, mutating generators, and Doxygen output.
- `make tests` and `make datagen` compile or link; they do not execute binaries. Record exact executable invocations separately.
- When HDF5 or source selection may change, inspect `configure`, `scripts/parser`, relevant recursive manifests, generated `SRC`, `TEXES`, `DGEXES`, and `IHDS`, and applicable default and `--disable-hdf5` builds. Do not assume all tabulated code disappears.
- Distinguish generator build, generator execution, fixture production, replay, perturbation comparison, checked return status, and assertion strength. A successful process with weak or discarded checks is not numerical proof.
- Run `.github/run_tests.sh` only in a disposable checkout: downloads, broad cleanup globs, and early-failure leftovers can affect pre-existing files.
- While `wiki/contradictions.md` records current `generate_makefile.sh` output as broken, inspect it only in a disposable tree and never run `make` from that output.
- Validate Doxygen with a unique temporary `OUTPUT_DIRECTORY` per `wiki/generated-boundaries.md`; do not place generated output beside `docs/raw/`.
- For Flux_Source or NRPyLeakage generated/derived C, establish the supported generator route first. Checked-in C and public headers own current behavior. Never invent commands or claim reproducibility when repo evidence records drift or unknown entry points.
- Workflow presence proves selection, not historical execution. Account for path filters, event types, compiler, OS, and HDF5 variants.
- Upstream `Unit_Tests/`, ET_Legacy tests, and static GRHayLib parity do not prove Cactus build, schedule, parameter, lifecycle, or runtime behavior. State when a real Cactus or Einstein Toolkit environment was unavailable.
- For public API changes, trace declarations, definitions, manifests, installed exposure, mode guards, dispatch, direct consumer compilation/linking, docs, tests, and downstream impact as applicable.
- For KB changes, keep repo-relative links, ground claims in exact files, keep routers compact, and run applicable checks from `wiki/lint/CHECKS.md`.

## Independent production

In `design`, each seat returns one complete design covering targets, interfaces, assumptions, existing components to reuse, risks, compatibility, generated ownership, and validation; end with `PROPOSAL COMPLETE`. In `tri-build`, each seat creates one complete candidate only in its assigned area and reports inspected inputs, rationale, reused machinery, validation, and limitations; end with `DRAFT COMPLETE`.

Freeze all results before comparison. Select one coherent proposal or candidate, or synthesize only compatible evidence-supported elements. Never average incompatible behavior or include an element merely to represent every seat. If none satisfies the brief, return `NO ACCEPTABLE CANDIDATE`. Retire proposal/build seats; use a fresh triad for candidate review.

## Review and decision

After required pre-review validation passes, each reviewer independently reports inventory and assumptions, role-specific evidence, every material finding or explicit none, and final exact line `DECISION: ACCEPT` or `DECISION: BLOCK`.

Each finding gives exact location, violated criterion or authority, evidence, consequence, and smallest sufficient correction or missing evidence. Unsupported preference, enlarged scope, conditional approval, silence, timeout, or stale-candidate approval does not pass.

Root classifies findings as `CONTRACT DEFECT`, `BLOCKER`, `EVIDENCE GAP`, `NONBLOCKING TRADEOFF`, or `NOISE`. After cycle three only, root may classify a residual as `MINOR FINALIZATION DEFECT`. Evidence decides classification, not vote count. A substantive block is not agent failure and may not be resampled. All three seats must return valid `ACCEPT` for the same frozen candidate before normal installation.

After cycle three there is no fourth review. One terminal correction batch is allowed only when every residual is a `MINOR FINALIZATION DEFECT`: local, deterministic, behavior-neutral, prescribed by existing authority, evidence-complete, and conclusively covered by existing checks. It must not change physics, algorithms, numerical behavior, public interfaces, manifests, dispatch, dependencies, compatibility, generated ownership, test meaning, KB claims or evidence relationships, or downstream contracts. Freeze and validate once. Report `LEAD-FINALIZED AFTER CYCLE CAP` and state that terminal bytes were not unanimously reviewed. Otherwise report `NO ACCEPTABLE CANDIDATE`.

## Install and report

Immediately before live writes, confirm authorization, unchanged baselines and decision-relevant inputs, unchanged candidate identity, and unanimous acceptance or valid terminal finalization. Apply only authorized operations. Verify installed bytes, object types, permissions, symlinks, required absences, and absence of validator artifacts; rerun affected installed-tree checks. If installation is partial or verification fails, stop and report exact live state; use only an already-authorized recovery path.

Report mode, scope/cycle IDs, wave/cycle/correction counts, recovery or transient-retry use, targets changed, validation and install verification, unavailable checks and untested configurations, reused GRHayL machinery, nonblocking tradeoffs, seat decisions bound to candidate identity, isolation and temporary-artifact status, exact blocker when nothing was installed, and outcome: `UNANIMOUSLY ACCEPTED`, `LEAD-FINALIZED AFTER CYCLE CAP`, or `NO ACCEPTABLE CANDIDATE`.
