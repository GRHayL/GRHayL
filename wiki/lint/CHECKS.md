# Manual KB Checks

These checks help maintain the GRHayL agent KB. Run them from repo root. They
are shell examples only, not a maintained script.

## Broken Repo-Relative Markdown Links

List Markdown links in KB pages:

```bash
rg -n '\[[^]]+\]\([^)]+\)' AGENTS.md wiki
```

Extract likely repo-relative Markdown link targets. Ignore external URLs,
anchors, and Doxygen `@ref` links in source docs:

```bash
rg -n -o '\[[^]]+\]\(([^)#][^)]*)\)' AGENTS.md wiki \
  | sed -E 's/.*\]\(([^)#]+).*/\1/' \
  | rg -v '^[a-z]+://|^mailto:' \
  | sort -u
```

Check each extracted target with `test -e`, resolving links relative to the
page that contains them:

```bash
rg -n -o '\[[^]]+\]\(([^)#][^)]*)\)' AGENTS.md wiki \
  | rg -v '\]\([a-z]+://|\]\(mailto:' \
  | while IFS=: read -r page line match; do
      target=$(printf '%s\n' "$match" | sed -E 's/.*\]\(([^)#]+).*/\1/')
      case "$target" in
        /*|../*) printf '%s:%s: review non-local target %s\n' "$page" "$line" "$target" ;;
        *) test -e "$(dirname "$page")/$target" || printf '%s:%s: missing %s\n' "$page" "$line" "$target" ;;
      esac
    done
```

During parallel KB work, do not fail solely because links to other agents'
pages are temporarily missing. Report them.

## Forbidden Policy Terms

Search for terms tied to banned source tracking. Hits in policy statements are
allowed only when they say the practice is forbidden.

```bash
rg -n '(sha[0-9]*|md5|hash(es|ing)?|checksum|digest|mtime|fingerprint|maintenance log|separate log)' AGENTS.md wiki || true
```

Inspect each hit:

```bash
rg -n -C 2 '(sha[0-9]*|md5|hash(es|ing)?|checksum|digest|mtime|fingerprint|maintenance log|separate log)' AGENTS.md wiki || true
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
find wiki -name '*.md' -maxdepth 4 -print | sort
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

- every KB page should be reachable from `wiki/index.md`, `wiki/catalog.md`,
  or a task-specific router page
- pages owned by other agents may be linked before they exist
- do not create placeholder files for pages owned by other agents

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
