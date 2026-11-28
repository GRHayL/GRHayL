# Knowledge Base Layout And Rules

> Page layout, source provenance, router bounds, and checker rules for the GRHayL Markdown KB. · Status: confirmed

## Summary

The GRHayL knowledge base separates authoritative repository sources (`GRHayL/`, `implementations/GRHayLib/`, `Unit_Tests/`, `docs/raw/`, `tools/compose/`, `scripts/`, `configure`, `generate_makefile.sh`, `Doxyfile`, `DoxygenLayout.xml`, `README.md`, `.github/`, and `raw/SOURCES.md`), compiled `wiki/` pages, and governance rules (`AGENTS.md`, `SCHEMA.md`, `workflows.md`). Navigation starts at `AGENTS.md`, branches through pure `index.md` routers (`2` to `14` children), and records each durable technical fact once in its defining leaf.

## Volatile Information Policy

Source-manifest tables in [raw/SOURCES.md](../raw/SOURCES.md) use `Source | Provenance | Status | Ingest` for aggregate rows, `Source | Status` for exact local file rows (added only when no [source-map.md](source-map.md) `Glob pattern` covers the path), and `Source | Tier | Provenance | Status | Ingest | Notes` (`<= 80` B) for external URL rows, where `Status` is `living`, `frozen`, or `frozen-historical`, `Ingest` is `registered`, `partial`, or `ingested`, and `Tier` is `external-spec` or `background`. Leaf `Status:` values (`confirmed`, `provisional`, `contested`, `stale`) and `raw/SOURCES.md` `Status`/`Ingest` states record current structural source-alignment state, never dated audit verdicts or test-run outcomes.

Authored KB files (`AGENTS.md`, `wiki/**/*.md`, and `raw/SOURCES.md`) forbid volatile snapshots: no dates, times, timestamps, `mtime` values, commit or file hashes/digests, inventory/file/page/source/test/job/function/line counts, line numbers or line ranges in citations, audit/reconciliation/resolution dates, environment or tool-version tuples, recorded run or review results, `"as of"` statements, or separate maintenance logs (Git history records changes). Stable scientific constants, array extents, stencil widths, tensor component counts, documented tolerances, publication years, versioned specification URLs, arXiv/DOI identifiers, durable code measurements (accuracy/convergence comparisons, memory footprints, and benchmark timings), and reviewed technical hash facts remain allowed.

Universal layered availability and execution rule: declaration in a header, C definition, `make.code.defn`/`configure` selection (default and `--disable-hdf5`), `make install` header copy, runtime EOS dispatch, `make tests`/`make datagen` compilation, binary execution, CI workflow selection, and external Cactus builds are distinct evidence layers. Do not repeat generic layered-availability disclaimers across individual leaves; document only module-specific unbuilt, stubbed, unwired, or untested paths. For NRPy-generated `Flux_Source` C files and ancestral `NRPyLeakage` code, checked-in C and headers govern compiled behavior.

## Support Pages

Every table row across `catalog.md`, `tests/unit-test-catalog.md`, `glossary.md`, `source-map.md`, `contradictions.md`, and `raw/SOURCES.md` must stay `<= 450` B:

- `wiki/catalog.md` and `wiki/tests/unit-test-catalog.md`: `<= 240` B/row; `catalog.md` uses the 3-column leaf index (`Page | Scope | Query terms`) listing each normal leaf once with non-duplicate symbol and topic search terms (no `Type` or `Status` columns; update it only when adding, moving, renaming, or deleting a leaf).
- `wiki/source-map.md`: aggregate scope table (`Source or aggregate | Glob pattern | Authority tier | Uncovered scope`, `<= 450` B/row) covering only the multi-file aggregates from `raw/SOURCES.md`, which alone records their provenance, `Status`, and `Ingest`; the `Glob pattern` column alone defines their path membership. Single-file rows and manual dependent-page lists are prohibited; compute reverse citations dynamically with `rg -l -F '<path>' wiki/`.
- `wiki/contradictions.md`: active conflict register (`ID | Claim status | Affected pages | Summary`, `<= 220` B/row) paired with self-contained `### CONTR-XXXX` subsections.
- `wiki/glossary.md`: canonical `Term | Meaning` table (`<= 450` B/row, `Meaning <= 180` B) where every `Meaning` cell links its defining leaf or marks external/background context, carries no `## Sources` section, and uses typed `## See Also` bullets.
- `wiki/workflows.md` and `wiki/computational-physics-language.md`: KB procedures, deterministic and review checks, and direct computational-physics prose rules.

## Page Layout

Normal leaves stay `<= 20,000` B and use this exact structure: `# Title`, two-line blockquote (`> One-line scope. · Status: <confirmed|provisional|contested|stale>` and `> Up: [Parent](...)`), `## Summary`, `## Detail`, `## Sources`, and `## See Also`.

- **Subsection partitioning and line bounds**: Whenever `## Detail` exceeds `4,000` B, partition it with descriptive `###` headings naming concrete functions, structs, parameters, or equations so no contiguous `## Detail` or `###` chunk exceeds `6,000` B (measured on raw UTF-8 text excluding HTML comments). Keep every `wiki/` line `<= 1,000` B, avoid double blank lines (`\n\n\n`), and never restate table cells in adjacent prose paragraphs.
- **Parent and neighbor links**: Header `> Up:` is the sole parent-router link (`../AGENTS.md` for leaves directly in `wiki/`, otherwise `index.md`). `## See Also` holds `1` to `6` non-parent neighbor `- ` bullets (targeting `<= 3` high-signal links, `1` unique link per bullet, no self-links), never duplicating `> Up:`, each starting with `Depends on:`, `Implements:`, `Validated by:`, `Example:`, `Contrasts with:`, or `See also:`.
- **Page status**: Use `contested` when an unresolved source conflict alters `Summary`, public interface, or multiple core `Detail` claims; `stale` when living-source drift leaves principal claims unreconciled; `provisional` when synthesis is incomplete without known conflict; and `confirmed` when registered sources directly support the page (allowing at most one bounded, noncentral `Detail` defect marked inline with its `CONTR-XXXX` link).

## Citation Rules

Source authority precedence is `primary-code` > `primary-test` > `generated-evidence` > `primary-doc` > `external-spec` > `background`; exceptions require explicit justification in the defining leaf or [contradictions.md](contradictions.md).

- **Unique `## Sources` entries**: List each cited repository file or external URL at most once as a clickable repo-relative or full official Markdown link in `## Sources`, grouping multiple functions, structs, macros, or headings on that bullet. `confirmed` pages must not cite `frozen-historical` sources unless the bullet explicitly notes historical provenance.
- **Doxygen coexistence**: Link `docs/raw/*.dox` and `docs/raw/derivation.md` as `primary-doc` instead of copying `@defgroup`, `@ingroup`, `@brief`, `@details`, `\f[`, or full API tables into `wiki/`. Keep `docs/raw/Reconstruction.dox` read-only during KB updates.
- **Single Defining-Page Rule**: State and evidence each repository fact once in its defining leaf. Neighboring leaves link to that defining leaf in one sentence rather than copying claims or `Evidence` lines.

### Claim And Evidence Rules

Legacy multi-line `Claim evidence:` blocks are prohibited. Add a single-line evidence tag (`<= 300` B) immediately after new or materially revised high-risk claims (public interfaces, mathematical/numerical guarantees, CLI options, generated-file boundaries, CI job behavior, authority/contradiction decisions, and easily misread source passages):

```text
Evidence (<role>): `<subpath.c::symbol>`; corroboration: `<source>`
```

Allowed roles are `descriptive behavior`, `normative rule`, `scientific/interface requirement`, `CI behavior`, and `generated evidence`. Cite primary source files (outside `wiki/` and `raw/SOURCES.md`) in `## Detail` via `` `subpath.c::symbol` `` tags matching `## Sources`, using the bare basename when no other tracked source file shares it and otherwise a repository path that names one file (for example `GRHayL/Flux_Source/make.code.defn`, `tools/compose/README.md`, or root `README.md` when no other `README.md` is in `## Sources`), omitting `corroboration:` when none exists. Inside tables, attach compact ``(Ev: `file.c::symbol`)`` tags directly in cells instead of repeating table rows in prose. Cite C struct members with dot or double-colon notation (for example `ghl.h::ghl_eos_parameters.rho_atm`), never pointer arrows (`->`). Cite binary HDF5 files (`.h5`) without `::symbol`.

## Router Rules

Pure directory routers (`wiki/**/index.md`) contain only a title, one-line blockquote, `Page | Go here when...` table, and parent/root links—never `## Detail`, `## Sources`, or standalone implementation facts. `AGENTS.md` (`<= 4,800` B) is the sole root router and pairs compact agent, language, and scope rules with the root `Task | Read first` router table. Every `wiki/**/index.md` router links `2` to `14` immediate children (leaves or child `index.md` routers); single-child padding directories are prohibited. `AGENTS.md` may link directly to `wiki/*.md` support pages, branch `wiki/*/index.md` routers, and interior `wiki/**/index.md` sub-routers, while every normal leaf and `index.md` sub-router must also be linked from its immediate parent router.

## Contradiction Rules

Active rows in [contradictions.md](contradictions.md) `## Register` use the 4-column header:

`ID | Claim status | Affected pages | Summary`

Rows stay `<= 220` B with `ID` in `CONTR-0001` form and `Claim status` in `{"contested", "stale"}`. Each active ID has a matching `### CONTR-XXXX` subsection recording competing sources, authority decision, `Evidence (<role>): ...` (`<= 300` B), recheck trigger, and `Resolution check:`. Every active affected leaf carries literal `Claim status: contested; contradiction: CONTR-0001.` (or `stale`) outside code fences, paired with a separate link to `contradictions.md#contr-0001`. When the resolution check passes and affected leaves are reconciled, remove all `#contr-*` links, delete the register row and `### CONTR-*` subsection in the same commit, and append the ID to `Retired IDs:` in `contradictions.md` (retired IDs are never reused).

## Checker And Frozen Scope

`PYTHONDONTWRITEBYTECODE=1 python3 -B docs/kb_lint.py` checks only that local links and heading anchors resolve in `AGENTS.md`, `wiki/**/*.md`, and `raw/SOURCES.md` (no wikilinks), that every wiki page is reachable from `AGENTS.md`, and that every backticked token in an `Evidence (...)` line or `(Ev: ...)` span names exactly one git-tracked source file present on disk outside `wiki/` and `raw/SOURCES.md` (`## Sources` links break ties first, then exact repository path) whose cited `::` symbols appear in that text file (as whole identifiers in code files); stage new source files with `git add` before citing them.

Keep KB Markdown within the plain subset parsed by `docs/kb_lint.py`: use top-level fenced code blocks (not fences inside block quotes or list items), inline `[text](path#anchor)` links without link titles, and headings in plain words and inline code without emphasis or links. `docs/kb_lint.py` follows its top-of-file frozen design rules (three regex checks only, no doctests or self-tests); every other rule on this page and in [Language For Computational Physicists](computational-physics-language.md) is checked in review.

### Root Planning And Task Files

At repository root only, commissioned planning or task Markdown files whose lowercase basename starts with `plan`, `task`, or `tasks` (or matches an explicitly requested root plan filename, with no segment matching `audit`, `scan`, `report`, `log`, `latest`, `token`, `rank`, or `ranking` in singular, plural, or numbered form) are exempt from KB routing checks and must never be filed, staged, moved, or deleted as KB content unless the user directs.

## Not In The KB

Generated outputs and scratch files (`Makefile`, `build/`, `test/`, `lib/`, `docs/html/`, `docs/latex/`, `nrpy/`, `*.bin`, `*.asc`, downloaded `*.h5`, `*.o`, caches, scratch logs, token/audit/scan/latest reports, and uncommissioned planning or ranking files) remain excluded from `wiki/` and `raw/`.

## Sources

- [AGENTS.md](../AGENTS.md) - root router, computational-physics prose mandate, and protected scope
- [docs/kb_lint.py](../docs/kb_lint.py) - link, reachability, and Evidence-citation checks

## See Also

- Depends on: [Language For Computational Physicists](computational-physics-language.md)
- See also: [Workflows](workflows.md)
- See also: [Glossary](glossary.md)
