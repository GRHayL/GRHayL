# Workflows

> Procedures for querying, ingesting, reconciling, reproducing, and checking the GRHayL KB. · Status: confirmed

## Summary

Start at [AGENTS.md](../AGENTS.md), follow branch and sub-branch `index.md` routers directly, and read only the needed leaf or `###` subsection. Record each durable fact once in its defining leaf, compute reverse file dependencies dynamically with `rg -l -F`, keep builds and mutating scripts inside disposable trees, and verify edits with `tools/kb_lint.py` and [Language For Computational Physicists](computational-physics-language.md).

## Search Order

1. Read [AGENTS.md](../AGENTS.md) and jump to the target branch `index.md` router (or `wiki/gems/atmosphere.md`).
2. For ambiguous or cross-branch queries, run `rg -n -i -F '<term>' wiki/catalog.md` and read only matching rows.
3. Read compact leaves (`< 8,000` B) directly in one read. On leaves `>= 8,000` B, run `rg -n '^### ' <page>` first and read only the needed `###` range. Follow a `## See Also` link only when its prefix (`Depends on:`, `Implements:`, `Validated by:`, `Example:`, `Contrasts with:`, `See also:`) identifies a missing prerequisite or check.
4. Search support tables selectively without reading whole files:
   - Reverse file citations: `rg -l -F '<basename>' wiki/`
   - Source or aggregate rows: ``rg -n -F '| `<key>` |' raw/SOURCES.md`` or ``rg -n -F '| `<aggregate>` |' wiki/source-map.md``
   - Contradiction entries: `rg -n -F '| <CONTR-ID> |' wiki/contradictions.md` plus the matching `### CONTR-XXXX` subsection
   - Canonical term spellings: `rg -n -i -F '<term>' wiki/glossary.md`
5. Reopen primary files in `GRHayL/`, `implementations/GRHayLib/`, `Unit_Tests/`, `docs/raw/`, `tools/compose/`, `configure`, or `.github/` directly when verifying or updating code-backed claims.

## Ingest

1. Confirm source coverage via an aggregate `Glob pattern` in [source-map.md](source-map.md) (registered by key in [raw/SOURCES.md](../raw/SOURCES.md)), an exact local `raw/SOURCES.md` row (`Source | Status`), or an external row (`Source | Tier | Provenance | Status | Ingest | Notes`, `Notes <= 80` B). Never store hashes, dates, timestamps, `mtime` values, inventory counts, line numbers, or run results.
2. Select the single defining leaf for the topic under the appropriate branch router.
3. Synthesize durable facts in direct computational-physics prose, add `<= 300` B `Evidence (<role>):` lines using unambiguous `` `subpath.c::symbol` `` tags ([KB Rules](SCHEMA.md#claim-and-evidence-rules)) for high-risk claims (C struct members as `struct.field`, binary `.h5` without `::symbol`), and list each cited file or URL once in `## Sources`.
4. Partition `## Detail` with descriptive `###` headings whenever `## Detail` exceeds `4,000` B (`<= 6,000` B per chunk, `<= 20,000` B per leaf) and keep `1` to `6` relation-prefixed non-parent `## See Also` bullets.
5. Update the immediate parent `index.md` and [catalog.md](catalog.md) only when adding, moving, renaming, or deleting a leaf; when an aggregate's path membership changes, edit only its `Glob pattern` in [source-map.md](source-map.md), and revise its `raw/SOURCES.md` provenance only when the described content changes.

## Query Filing

Do not file one-off answers. File durable syntheses only when they resolve recurring cross-module questions, record non-obvious numerical or build constraints, or fix missing router paths:

1. Single-topic updates go into the existing defining leaf (splitting into two focused leaves if the leaf would exceed `20,000` B).
2. Narrow new topics add a leaf under the branch router (`2` to `14` children per `index.md`).
3. Cross-branch comparisons with one leaf sit directly under `wiki/` (`> Up: [GRHayL Knowledge Base](../AGENTS.md)`); create a dedicated directory router only when at least two related leaves exist.

## Living-Source Change

1. Inspect the changed code, header, manifest, Doxygen, or CI diff and verify source coverage (Ingest step 1) plus [Change Impact](architecture/change-impact.md).
2. Locate citing wiki pages dynamically with `rg -l -F '<changed_basename>' wiki/`: prose, tables, and Evidence tags often name files by basename or a partial path, and a basename search also matches full paths; for a file named inside a `{a,b}` group, also search the text through the brace (for example `compute_{`).
3. Run `rg -l -F '<changed_symbol>' wiki/` across other wiki pages only when a public function signature, enum/struct layout, header guard, `make.code.defn` manifest entry, or mathematical convention changed—never expand 1-hop across generic [glossary.md](glossary.md) terms.
4. Reconcile the defining leaf; replace duplicate prose in neighboring leaves with a one-sentence link. Mark a leaf `stale` with a `CONTR-XXXX` entry only when immediate reconciliation is blocked.
5. **Downstream GRHayLib coordination:** Before editing `implementations/GRHayLib/` or public `GRHayL/include/` structs/enums read by `GRHayLib`, check [GRHayLib Runtime Parameters](implementations/grhaylib-runtime-parameters.md) and [GRHayLib Verification](implementations/grhaylib-verification.md). Static parity in `implementations/GRHayLib/` does not prove Cactus build, schedule, or runtime behavior.

## External Source And Contradiction Updates

Classify external URLs in [raw/SOURCES.md](../raw/SOURCES.md) as `background` unless a durable KB claim depends on a version-pinned paper, tag, or format specification (`external-spec`). Repository code, tests, generated evidence, and repository docs take precedence per [KB Rules](SCHEMA.md#citation-rules).

When sources conflict or living-source drift cannot be reconciled in the same change, add a 4-column row (`<= 220` B) and `### CONTR-XXXX` subsection in [contradictions.md](contradictions.md) and place literal `Claim status: contested; contradiction: CONTR-0001.` (or `stale`) plus a `#contr-0001` link on every affected leaf. When resolved and verified by its resolution check, reconcile affected leaves, remove all `#contr-*` links, delete the register row and subsection, and append the ID to `Retired IDs:`.

## Safe Reproduction

In the shared repository checkout, run only side-effect-free scoped checks (`PYTHONDONTWRITEBYTECODE=1 python3 -B tools/kb_lint.py`, `git diff HEAD --check`, `./configure --help`, `./scripts/parser`, `sh -n`/`bash -n`, read-only YAML parsing). Run `./configure`, `make`, `make tests`, `make datagen`, `.github/run_tests.sh`, `generate_makefile.sh`, CompOSE converters, and Doxygen builds inside a unique disposable copy with redirected output paths (`OUTPUT_DIRECTORY` for Doxygen, never beside `docs/raw/`). Inspect scripts for downloads, cleanup globs, and failure leftovers before execution; never run `make` on broken `generate_makefile.sh` output, never edit `docs/raw/Reconstruction.dox` during KB work, and never add or enable tests or CI test executions without explicit user permission.

## Evidence-First Debugging

Follow the runtime failure inspection rule in [AGENTS.md](../AGENTS.md#agent-working-rules): reproduce or explain why not, inspect actual build outputs, generated `Makefile` variables (`SRC`, `TEXES`, `DGEXES`, `IHDS`), or runtime quantities before proposing a cause or fix, and record only durable lessons in `wiki/`.

## Deterministic And Review Checks

Run before completing KB edits:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -B tools/kb_lint.py
git diff HEAD --check
git diff --cached --check
```

`git diff HEAD --check` checks the working tree against `HEAD` and `git diff --cached --check` checks the staged copy; run `git add -N` on new pages first so the working-tree check sees them.

`tools/kb_lint.py` must exit `0` and print `KB lint passed.` In a shared checkout, if the failures come only from another author's uncommitted changes, run the check in a disposable clone that carries only your changes, which must pass, and report those failures without editing that author's files.

Optional review-time read-only scans (not part of `tools/kb_lint.py`):

```bash
# Forbidden volatile metadata, hashes, line citations, dates, and run results:
rg -n -i '(checksum|digest|mtime|fingerprint|maintenance log|\bas of\b|\.(c|h|py|sh|yml|md|dox|ccl):[0-9]+)' AGENTS.md wiki raw/SOURCES.md || true
# Doxygen markup duplication check:
rg -n 'defgroup|@ingroup|@brief|@details|\\f\[' wiki || true
# Absolute workspace path check (outside operational delivery-path clause in AGENTS.md):
rg -n '(^|[^[:alnum:]_.])/work(/|`|$)' wiki raw/SOURCES.md || true
```

Semantic review checklist:
- `git status --short` shows only intended changes, with no generated `Makefile`, `build/`, `test/`, `lib/`, `docs/html/`, `docs/latex/`, cache, or scratch report files added.
- `AGENTS.md <= 4,800` B; every normal leaf `<= 20,000` B with `## Detail` chunks `<= 6,000` B; every `wiki/` line `<= 1,000` B; every router links `2`–`14` children with no prose facts; support tables obey `SCHEMA.md` row byte caps.
- Every new or revised passage follows [KB Rules](SCHEMA.md) and [Language For Computational Physicists](computational-physics-language.md) with zero vague metaphors, no prose restating adjacent table cells, and single-defining-page authority.

## Sources

- [AGENTS.md](../AGENTS.md) - root router, working rules, computational-physics prose mandate, and test/CI restriction
- [tools/kb_lint.py](../tools/kb_lint.py) - link, reachability, and Evidence-citation checks
- [README.md](../README.md) - build, install, and test workflow overview
- [configure](../configure) - build configuration options and HDF5 flags

## See Also

- Depends on: [KB Rules](SCHEMA.md)
- Depends on: [Language For Computational Physicists](computational-physics-language.md)
- See also: [Build And Configure](architecture/build-and-configure.md)
- See also: [Change Impact](architecture/change-impact.md)
- See also: [Tests](tests/index.md)
