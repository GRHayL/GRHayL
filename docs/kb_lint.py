#!/usr/bin/env python3
"""
Minimalist link, reachability, and Evidence checker for the GRHayL Markdown KB.

FROZEN RULES:
- Scope: single standard-library Python script with ONLY the three regex checks
  below. Only fix defects in these three checks; do not add checks, CLI
  options, features, or third-party dependencies.
- No meta-linting: NO DOCTESTS ALLOWED (`doctest` is forbidden), no self-tests,
  no unit tests for this script, and no static-analysis pragmas or wrappers.
- Plain Markdown only (no CommonMark parser): mask only top-level fenced code
  blocks (0-3 spaces, ``` or ~~~, closing fence of >= opening length with the
  same marker), HTML comments (<!--...-->), and backtick inline code spans. Do
  not add parsers or regexes for blockquoted/list-indented fences, 4-space
  indented code, backslash-escaped backticks, reference links, link titles, or
  heading markup; write KB pages in the plain subset.
- Review-only rules stay out of code: page layout, byte limits, table columns,
  volatile-information prohibitions, and direct prose live in `wiki/SCHEMA.md`
  and `wiki/computational-physics-language.md`, not in this script.

Checks on `AGENTS.md`, `wiki/**/*.md`, and `raw/SOURCES.md`:
1. Local links & anchors: inline local links resolve inside the repository and
   `#anchor` matches a GitHub-style heading slug (`-1`, `-2` for duplicates);
   `[[wikilinks]]` fail.
2. Reachability: every `wiki/**/*.md` page is reachable from `AGENTS.md`.
3. Evidence citations: every backticked token in a wiki `Evidence (...)` line
   or `(Ev: ...)` span names one git-tracked file on disk outside `wiki/` and
   `raw/SOURCES.md` (ties break by `## Sources` links first, then exact repo
   path), and any `::symbol` exists in that text file (whole identifier in code
   files; substring in other text files; `.h5` files reject `::symbol`).

Adapted from NRPy's tools/kb_lint.py (Zachariah B. Etienne).
"""

import re
import subprocess
import sys
from collections import deque
from pathlib import Path
from urllib.parse import unquote, urlparse

ROOT = Path(__file__).resolve().parents[1]
WIKI = ROOT / "wiki"
AGENTS = ROOT / "AGENTS.md"
SOURCES = ROOT / "raw" / "SOURCES.md"

LINK_RE = re.compile(r"(?<!!)\[((?:[^\]\[]|\[[^\]\[]*\])+)\]\(([^)\n]+)\)")
INLINE_CODE_RE = re.compile(
    r"(?<!`)(?P<ticks>`+)(?!`)(?P<body>.*?)(?<!`)(?P=ticks)(?!`)", re.DOTALL
)
FENCE_OR_COMMENT_RE = re.compile(
    r"^[ ]{0,3}(`{3,})[^`\n]*\n.*?^[ ]{0,3}\1`*[ \t]*$"
    r"|^[ ]{0,3}(~{3,})(?!~)[^\n]*\n.*?^[ ]{0,3}\2~*[ \t]*$"
    r"|<!--.*?-->",
    re.DOTALL | re.MULTILINE,
)
CODE_SUFFIXES = (".c", ".h", ".cc", ".cpp", ".hpp", ".py")


def _fail(failures, path, line, message):
    rel = path.resolve().relative_to(ROOT).as_posix()
    failures.append(f"{rel}: {message}" if line is None else f"{rel}:{line}: {message}")


def _blank(match):
    return re.sub(r"[^\n]", " ", match.group(0))


def _prose(path):
    return FENCE_OR_COMMENT_RE.sub(_blank, path.read_text(encoding="utf-8"))


def _heading_anchors(path):
    anchors, counts = set(), {}
    for m in re.finditer(r"^#{1,6}[ \t]+(.+?)[ \t]*$", _prose(path), re.MULTILINE):
        title = re.sub(r"\s+#+\s*$", "", m.group(1)).strip().lower()
        base = re.sub(r"\s", "-", re.sub(r"[^\w\s-]", "", title))
        if not base:
            continue
        slug = base
        while slug in anchors:
            counts[base] = counts.get(base, 0) + 1
            slug = f"{base}-{counts[base]}"
        anchors.add(slug)
    return anchors


def _check_links(files, failures):
    graph = {}
    for path in files:
        text = INLINE_CODE_RE.sub(_blank, _prose(path))
        graph[path] = []
        for m in re.finditer(r"\[\[[^\]]+\]\]", text):
            _fail(failures, path, text.count("\n", 0, m.start()) + 1, "wikilink")
        for m in LINK_RE.finditer(text):
            line = text.count("\n", 0, m.start()) + 1
            target = m.group(2).strip()
            parsed = urlparse(target)
            if parsed.scheme in {"http", "https", "mailto", "tel"} or parsed.netloc:
                continue
            anchor = unquote(parsed.fragment)
            if not parsed.path:
                resolved = path
            elif parsed.path.startswith("/"):
                resolved = (ROOT / unquote(parsed.path).lstrip("/")).resolve()
            else:
                resolved = (path.parent / unquote(parsed.path)).resolve()
            if ROOT not in resolved.parents and resolved != ROOT:
                _fail(failures, path, line, f"link target escapes repository: {target}")
            elif not resolved.exists():
                _fail(failures, path, line, f"link target missing: {target}")
            elif anchor and (resolved.suffix != ".md" or anchor not in _heading_anchors(resolved)):
                _fail(failures, path, line, f"link anchor missing: {target}")
            else:
                graph[path].append(resolved)
    return graph


def _check_reachability(pages, graph, failures):
    root_doc = AGENTS.resolve()
    reachable = set()
    queue = deque([root_doc])
    while queue:
        path = queue.popleft()
        if path not in reachable:
            reachable.add(path)
            queue.extend(t for t in graph.get(path, []) if t == root_doc or WIKI in t.parents)
    for page in pages:
        if page not in reachable:
            _fail(failures, page, None, "wiki page not reachable from AGENTS.md")


def _symbol_problem(symbol, source_rel):
    if source_rel.endswith(".h5"):
        return f"cannot check Evidence symbol {symbol!r} in binary file {source_rel}"
    segments = [seg.strip() for seg in symbol.split("::")]
    valid = all(segments)
    is_code = source_rel.endswith(CODE_SUFFIXES)
    if is_code:
        base = symbol.split("(", 1)[0].strip()
        valid = valid and bool(re.fullmatch(r"[A-Za-z_]\w*(?:(?:\.|::)[A-Za-z_]\w*)*", base))
        segments = re.split(r"::|\.", base)
    if not valid:
        return f"invalid Evidence symbol syntax {symbol!r} in {source_rel}"
    try:
        source_text = (ROOT / source_rel).read_text(encoding="utf-8")
    except UnicodeDecodeError:
        return f"cannot check Evidence symbol {symbol!r} in binary file {source_rel}"
    if is_code:
        missing = [s for s in segments if not re.search(rf"(?<!\w){re.escape(s)}(?!\w)", source_text)]
    else:
        missing = [s for s in segments if s not in source_text]
    return f"Evidence symbol {missing[0]!r} not found in {source_rel}" if missing else None


def _check_evidence(pages, failures):
    listing = subprocess.run(
        ["git", "ls-files", "-z"], cwd=ROOT, capture_output=True, text=True, check=False
    )
    if listing.returncode != 0:
        raise RuntimeError(f"kb_lint needs a git checkout: {listing.stderr.strip()}")
    tracked = [
        p
        for p in listing.stdout.split("\0")
        if p and not p.startswith("wiki/") and p != "raw/SOURCES.md" and (ROOT / p).is_file()
    ]
    for page in pages:
        text = _prose(page)
        src_block = re.search(r"^## Sources\s*$(.*?)(?=^## |\Z)", text, re.M | re.S)
        sources = {
            (page.parent / unquote(m.group(2).split("#", 1)[0].strip())).resolve()
            for m in LINK_RE.finditer(src_block.group(1) if src_block else "")
        }
        for line_no, line in enumerate(text.splitlines(), start=1):
            stripped = line.strip()
            if stripped.startswith("Evidence ("):
                bodies = [re.sub(r"^Evidence \([^)]+\):", "", stripped, count=1)]
            else:
                masked = INLINE_CODE_RE.sub(_blank, stripped)
                bodies, pos = [], 0
                for m in re.finditer(r"\(Ev:", masked):
                    if m.start() < pos:
                        continue
                    depth, end = 1, m.end()
                    while end < len(masked) and depth:
                        depth += (masked[end] == "(") - (masked[end] == ")")
                        end += 1
                    bodies.append(stripped[m.end() : end - (depth == 0)])
                    pos = end
            for body in bodies:
                for m in INLINE_CODE_RE.finditer(body):
                    head, sep, symbol = (part.strip() for part in m.group("body").partition("::"))
                    matches = [t for t in tracked if t == head or t.endswith("/" + head)]
                    matches = (
                        [t for t in matches if ROOT / t in sources]
                        or [t for t in matches if t == head]
                        or matches
                    )
                    if len(matches) != 1:
                        msg = (
                            f"ambiguous Evidence file {head!r} matches {', '.join(matches[:3])}"
                            if matches
                            else f"Evidence file {head!r} is not a git-tracked file on disk"
                        )
                    else:
                        msg = _symbol_problem(symbol, matches[0]) if sep else None
                    if msg:
                        _fail(failures, page, line_no, msg)


def _main():
    pages = sorted(p.resolve() for p in WIKI.rglob("*.md"))
    files = sorted({AGENTS.resolve(), SOURCES.resolve(), *pages})
    failures = []
    graph = _check_links(files, failures)
    _check_reachability(pages, graph, failures)
    _check_evidence(pages, failures)
    if failures:
        print("KB lint failed:")
        for item in sorted(failures):
            print(f"- {item}")
        return 1
    print("KB lint passed.")
    return 0


if __name__ == "__main__":
    sys.exit(_main())
