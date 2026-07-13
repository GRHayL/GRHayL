# Manual KB Checks

These checks help maintain the GRHayL agent KB. Run them from repo root. They
are shell examples only, not a maintained script.

## Broken Repo-Relative Markdown Links

This dependency-free check resolves `../` links relative to each page, rejects
absolute/escaping paths, and reports missing targets:

```bash
python3 - <<'PY'
from pathlib import Path
import re

root = Path.cwd().resolve()
pages = [Path("AGENTS.md"), *sorted(Path("wiki").rglob("*.md"))]
link = re.compile(r"\[[^]]+\]\(([^)]+)\)")
for page in pages:
    for line_number, line in enumerate(page.read_text().splitlines(), 1):
        for raw in link.findall(line):
            target = raw.split("#", 1)[0]
            if not target or re.match(r"^[A-Za-z][A-Za-z0-9+.-]*:", target):
                continue
            candidate = (page.parent / target).resolve()
            try:
                candidate.relative_to(root)
            except ValueError:
                print(f"{page}:{line_number}: escaping target: {raw}")
                continue
            if not candidate.exists():
                print(f"{page}:{line_number}: missing target: {raw}")
PY
```

This is a conservative Markdown regex, not a full GitHub-Flavored Markdown
parser. Manually review links containing escaped/nested parentheses, reference
definitions, HTML, or titles; do not turn uncertain syntax into a hard failure.

If KB pages are being authored concurrently, do not fail solely because
links to in-progress pages are temporarily missing. Report them.

## Forbidden Policy Terms

Search for terms tied to banned source tracking. Hits in policy statements are
allowed only when they say the practice is forbidden.

```bash
rg -n '(\bsha([0-9]+)?(sum)?\b|\bmd5(sum)?\b|\bhash(es|ing)?\b|checksum|digest|\bmtime\b|fingerprint|maintenance log|separate log)' AGENTS.md wiki || true
```

Inspect each hit:

```bash
rg -n -C 2 '(\bsha([0-9]+)?(sum)?\b|\bmd5(sum)?\b|\bhash(es|ing)?\b|checksum|digest|\bmtime\b|fingerprint|maintenance log|separate log)' AGENTS.md wiki || true
```

Policy must remain:

- no source-tracking hashes or hashing of sources
- no `mtime`
- no stored fingerprints
- no separate maintenance logs
- dates only when absolutely necessary, and retained dates use `MM-DD-YYYY`

Find date-like text for manual review:

```bash
rg -n '([0-9]{4}-[0-9]{2}-[0-9]{2}|[0-9]{2}[/-][0-9]{2}[/-][0-9]{4}|[A-Z][a-z]+ [0-9]{1,2}, [0-9]{4})' AGENTS.md wiki || true
```

## Page Contract Checks

Each KB page should state purpose, link ground truth, and avoid replacing
Doxygen or source.

```bash
find wiki -maxdepth 4 -name '*.md' -print | sort
```

For each page, verify:

- purpose appears near the top
- repo-relative source/doc/test paths support substantive claims
- `Ground Truth References` appears only if web facts were used
- external references, if any, are official full URLs
- no absolute workspace paths
- no Doxygen table copied wholesale when a route would suffice

Search for absolute workspace paths:

```bash
rg -n '(^|[^[:alnum:]_.])/work(/|`|$)' AGENTS.md wiki || true
```

## Orphan Page Checks

List KB pages and links to KB pages:

```bash
find wiki -name '*.md' -print | sort
rg -n '\]\(([^)]*wiki/|[^)]*\.md)' AGENTS.md wiki || true
```

Manual review:

- every KB page should be reachable from `wiki/index.md` or `wiki/catalog.md`
- links to planned-but-unwritten pages may exist temporarily; report them
- do not create placeholder files for pages you are not writing

## Doxygen Duplication Checks

KB pages should route to Doxygen source, not duplicate long Doxygen passages.
Review large pages and repeated Doxygen headings:

```bash
wc -l wiki/*.md wiki/lint/*.md 2>/dev/null
rg -n 'defgroup|@ingroup|@brief|@details|\\f\\[' wiki || true
```

If a page repeats equations, solver tables, or function docs from `docs/raw/`,
replace the copy with a link to the relevant Doxygen source:

```bash
find docs/raw -maxdepth 1 -type f | sort
```

## Build, Workflow, And Generated-Boundary Checks

Read-only option/parser checks:

```bash
./configure --help
./scripts/parser awk GRHayL/make.code.defn subdirs
./scripts/parser awk GRHayL/include/make.code.defn install_headers
find GRHayL -name make.code.defn -print | sort
find Unit_Tests -maxdepth 1 -name 'unit_test_*.c' -print | sort
find Unit_Tests/data_gen -maxdepth 1 -name 'unit_test_data_*.c' -print | sort
```

Compare help against `configure`'s build-type `case` manually. Current known
faults are help-only `nocflags`, parser-only `plain`, and differing production
flag strings; do not call either spelling supported without reopening source.

Check shell syntax without execution:

```bash
sh -n configure scripts/parser
bash -n generate_makefile.sh .github/run_tests.sh
```

Parse workflow YAML when PyYAML is available. `BaseLoader` preserves the key
`on` as text instead of applying YAML 1.1 boolean conversion:

```bash
python3 - <<'PY'
from pathlib import Path
import yaml

for path in sorted(Path(".github/workflows").glob("*.yml")):
    data = yaml.load(path.read_text(), Loader=yaml.BaseLoader)
    print(path, "events=", sorted(data["on"]), "jobs=", len(data["jobs"]))
PY
```

Workflow/job presence is selection evidence, not execution evidence. Compare
job commands with `.github/actions/*/action.yml`, `.github/run_tests.sh`, and
actual generated Makefile targets before claiming compile, run, or coverage.

Run Doxygen only with the unique temporary-output procedure in
[Generated Boundaries](../generated-boundaries.md). Record `doxygen --version`
and compare warnings only with the same version/configuration. Confirm no
generated workspace outputs remain:

```bash
find . -maxdepth 2 \( -name Makefile -o -name build -o -name test \
  -o -path './docs/html' -o -path './docs/latex' -o -path './lib' \) -print
```

`generate_makefile.sh` is currently broken. Reproduce it only in a unique
disposable checkout, inspect generated variables/targets, and never invoke
`make` when malformed targets, outside-`GRHayL` sources, or invalid include
roots are present.

## Source Authority Checks

When source paths change, review affected KB pages by dependency, not by stored
source fingerprints.

```bash
git status --short
git diff --name-only
```

Map changed paths to KB review targets:

```bash
git diff --name-only | while read -r path; do
  case "$path" in
    GRHayL/include/*) echo "$path -> wiki/public-api-map.md wiki/catalog.md" ;;
    GRHayL/Atmosphere/*) echo "$path -> wiki/gems/index.md wiki/catalog.md" ;;
    GRHayL/Con2Prim/*) echo "$path -> wiki/gems/index.md wiki/catalog.md" ;;
    GRHayL/EOS/*) echo "$path -> wiki/gems/index.md wiki/catalog.md wiki/build-and-ci.md" ;;
    GRHayL/Flux_Source/*) echo "$path -> wiki/physics/evolution-equation-map.md wiki/catalog.md" ;;
    GRHayL/Induction/*) echo "$path -> wiki/physics/evolution-equation-map.md wiki/catalog.md" ;;
    GRHayL/Neutrinos/*) echo "$path -> wiki/gems/index.md wiki/catalog.md" ;;
    GRHayL/Reconstruction/*) echo "$path -> wiki/gems/index.md wiki/catalog.md" ;;
    Unit_Tests/*) echo "$path -> wiki/test-map.md wiki/catalog.md" ;;
    docs/raw/*) echo "$path -> wiki/index.md wiki/catalog.md affected topic page" ;;
    .github/*|configure|generate_makefile.sh) echo "$path -> wiki/build-and-ci.md wiki/generated-boundaries.md" ;;
  esac
done
```

Before editing any KB page, re-open the named source, header, Doxygen source,
test, or CI file and make the smallest source-backed update.
