#!/usr/bin/env python3
"""
Link, reachability, and citation checks for the GRHayL Markdown knowledge base.

DO NOT ADD NEW CHECKS OR FEATURES TO THIS FILE. Only fix defects in the three
existing checks below. KB structure, size, and prose rules are review guidance
in wiki/SCHEMA.md and wiki/computational-physics-language.md, not lint checks.

The checks scan AGENTS.md, wiki/**/*.md, and raw/SOURCES.md, ignoring fenced
code, HTML comments, and inline code:

1. Every local Markdown link resolves to an existing repository file and, when
   it has an anchor, to an existing heading; wikilinks are rejected.
2. Every wiki page is reachable from AGENTS.md through links among AGENTS.md
   and wiki pages.
3. Every backticked source file in a wiki ``Evidence (...)`` line or
   ``(Ev: ...)`` span names exactly one git-tracked file present on disk
   (``## Sources`` links break ties), and every cited ``file::symbol``
   target appears in that file.

Adapted from NRPy's tools/kb_lint.py (Zachariah B. Etienne).
"""

# Regexes are needed for Markdown link, heading, and code-span parsing.
import re
import subprocess
import sys
from collections import deque
from pathlib import Path
from typing import Deque, Dict, List, Match, Optional, Set, Tuple
from urllib.parse import unquote, urlparse

ROOT = Path(__file__).resolve().parents[1]
WIKI = ROOT / "wiki"
AGENTS = ROOT / "AGENTS.md"
SOURCES = ROOT / "raw" / "SOURCES.md"
LINK_RE = re.compile(r"(?<!!)\[([^\]\[]+)\]\(([^)\n]+)\)")
INLINE_CODE_RE = re.compile(
    r"(?<!`)(?P<ticks>`+)(?!`)(?P<body>.*?)(?<!`)(?P=ticks)(?!`)", re.DOTALL
)
FENCE_OR_HTML_COMMENT_RE = re.compile(
    r"^[ ]{0,3}`{3,}[^\n]*\n.*?^[ ]{0,3}`{3,}[ \t]*$"
    r"|^[ ]{0,3}~{3,}[^\n]*\n.*?^[ ]{0,3}~{3,}[ \t]*$"
    r"|<!--.*?-->",
    re.DOTALL | re.MULTILINE,
)
EVIDENCE_FILE_SUFFIX_RE = re.compile(
    r"(?:\.(?:c|h|cc|cpp|hpp|py|sh|awk|yml|yaml|ccl|dox|tex|bib|css|xml|md|txt|h5|defn|cfg|ini)$"
    r"|(?:^|/)(?:Makefile|Doxyfile|configure|parser|README|AUTHORS|LICENSE"
    r"|THIRD_PARTY_NOTICES|\.clang-format|\.editorconfig|\.gitattributes"
    r"|\.gitignore|con2prim-testdata-ref|et-legacy-testdata-ref)$)"
)
CODE_SUFFIXES = (".c", ".h", ".cc", ".cpp", ".hpp", ".py")
BINARY_SUFFIXES = (".h5",)


def _fail(failures: List[str], path: Path, line: Optional[int], message: str) -> None:
    """
    Append a formatted lint failure.

    :param failures: Mutable failure list.
    :param path: File associated with the failure.
    :param line: Optional one-based line number.
    :param message: Failure detail.
    """
    rel = path.resolve().relative_to(ROOT).as_posix()
    failures.append(f"{rel}: {message}" if line is None else f"{rel}:{line}: {message}")


def _blank(match: Match[str]) -> str:
    """
    Replace a match with spaces, keeping newlines so offsets and line numbers hold.

    :param match: Regex match to blank.
    :return: Same-length replacement text.
    """
    return re.sub(r"[^\n]", " ", match.group(0))


def _prose(path: Path) -> str:
    """
    Return file text with fenced code and HTML comments blanked.

    :param path: Markdown file to read.
    :return: Masked text with the same character offsets as the file.
    """
    return FENCE_OR_HTML_COMMENT_RE.sub(_blank, path.read_text(encoding="utf-8"))


def _without_inline_code(text: str) -> str:
    """
    Blank Markdown inline-code spans while keeping character offsets.

    :param text: Masked Markdown text.
    :return: Text with inline-code spans blanked.
    """
    return INLINE_CODE_RE.sub(_blank, text)


def _heading_anchors(path: Path) -> Set[str]:
    """
    Return GitHub-style anchor slugs for the Markdown headings in a file.

    :param path: Markdown file to scan.
    :return: Anchor slugs.
    """
    anchors: Set[str] = set()
    for match in re.finditer(r"^#{1,6}[ \t]+(.+?)[ \t]*$", _prose(path), re.MULTILINE):
        title = re.sub(r"\s+#+\s*$", "", match.group(1)).strip().lower()
        anchors.add(re.sub(r"\s+", "-", re.sub(r"[^\w\s-]", "", title)))
    return anchors


def _check_links(files: List[Path], failures: List[str]) -> Dict[Path, List[Path]]:
    """
    Check local links, anchors, and wikilinks; return the local link graph.

    :param files: Markdown files to scan.
    :param failures: Mutable failure list.
    :return: Map from each file to the existing local files it links to.
    """
    graph: Dict[Path, List[Path]] = {}
    for path in files:
        text = _without_inline_code(_prose(path))
        graph[path] = []
        for match in re.finditer(r"\[\[[^\]]+\]\]", text):
            _fail(failures, path, text.count("\n", 0, match.start()) + 1, "wikilink")
        for match in LINK_RE.finditer(text):
            line = text.count("\n", 0, match.start()) + 1
            target = match.group(2).strip()
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
            elif anchor and (
                resolved.suffix != ".md"
                or re.sub(r"\s+", "-", re.sub(r"[^\w\s-]", "", anchor.lower()))
                not in _heading_anchors(resolved)
            ):
                _fail(failures, path, line, f"link anchor missing: {target}")
            else:
                graph[path].append(resolved)
    return graph


def _check_reachability(
    pages: List[Path], graph: Dict[Path, List[Path]], failures: List[str]
) -> None:
    """
    Check that every wiki page is reachable from ``AGENTS.md`` through local links.

    :param pages: Wiki pages.
    :param graph: Local link graph from ``_check_links``.
    :param failures: Mutable failure list.
    """
    reachable: Set[Path] = set()
    queue: Deque[Path] = deque([AGENTS.resolve()])
    while queue:
        path = queue.popleft()
        if path not in reachable:
            reachable.add(path)
            queue.extend(
                t
                for t in graph.get(path, [])
                if t == AGENTS.resolve() or WIKI in t.parents
            )
    for page in pages:
        if page not in reachable:
            _fail(failures, page, None, "wiki page not reachable from AGENTS.md")


def _evidence_tokens(span: str) -> List[Tuple[str, bool]]:
    """
    Return ``(token, inside_parentheses)`` for each backticked token in an Evidence span.

    :param span: Evidence line or ``(Ev: ...)`` substring.
    :return: Tokens in order, flagged when nested inside parentheses.
    """
    header = re.match(r"^(?:Evidence \([^)]+\):|\(Ev:)", span)
    body = span[header.end() :] if header else span
    tokens: List[Tuple[str, bool]] = []
    depth = 0
    idx = 0
    while idx < len(body):
        if body[idx] == "`":
            code_match = INLINE_CODE_RE.match(body, idx)
            if code_match:
                raw = code_match.group("body")
                if raw.startswith(" ") and raw.endswith(" ") and raw.strip(" "):
                    raw = raw[1:-1]
                tokens.append((raw, depth > 0))
                idx = code_match.end() - 1
            else:
                idx += len(body[idx:]) - len(body[idx:].lstrip("`")) - 1
        elif body[idx] == "(":
            depth += 1
        elif body[idx] == ")":
            if depth > 0:
                depth -= 1
            elif span.startswith("(Ev:"):
                break
        idx += 1
    return tokens


def _symbol_problem(symbol: str, source_rel: str) -> Optional[str]:
    """
    Return why a cited ``::`` or ``struct.field`` symbol fails, if it does.

    :param symbol: Symbol text after ``::`` or the whole parenthesized token.
    :param source_rel: Repository-relative path of the cited tracked source file.
    :return: Failure message, or None when every cited segment is in the file.
    """
    if source_rel.endswith(BINARY_SUFFIXES):
        return f"cannot check Evidence symbol {symbol!r} in binary file {source_rel}"
    segments = [seg.strip() for seg in symbol.split("::")]
    valid = all(segments)
    if source_rel.endswith(CODE_SUFFIXES):
        base = symbol.split("(", 1)[0].strip()
        valid = valid and bool(
            re.fullmatch(r"[A-Za-z_]\w*(?:(?:\.|::)[A-Za-z_]\w*)*", base)
        )
        segments = re.split(r"::|\.", base)
    if not valid:
        return f"invalid Evidence symbol syntax {symbol!r} in {source_rel}"
    try:
        source_text = (ROOT / source_rel).read_text(encoding="utf-8")
    except UnicodeDecodeError:
        return f"cannot check Evidence symbol {symbol!r} in binary file {source_rel}"
    missing = [seg for seg in segments if seg not in source_text]
    return (
        f"Evidence symbol {missing[0]!r} not found in {source_rel}" if missing else None
    )


def _check_evidence(pages: List[Path], failures: List[str]) -> None:
    """
    Check that Evidence citations name one git-tracked file and existing symbols.

    :param pages: Wiki pages.
    :param failures: Mutable failure list.
    :raises RuntimeError: If ``git ls-files`` fails.
    """
    listing = subprocess.run(
        ["git", "ls-files", "-z"], cwd=ROOT, capture_output=True, text=True, check=False
    )
    if listing.returncode != 0:
        raise RuntimeError(f"kb_lint needs a git checkout: {listing.stderr.strip()}")
    tracked = [
        p
        for p in listing.stdout.split("\0")
        if p
        and not p.startswith("wiki/")
        and p != "raw/SOURCES.md"
        and (ROOT / p).is_file()
    ]
    for page in pages:
        text = _prose(page)
        sources_match = re.search(r"^## Sources\s*$(.*?)(?=^## |\Z)", text, re.M | re.S)
        sources = {
            (page.parent / unquote(m.group(2).split("#", 1)[0].strip())).resolve()
            for m in LINK_RE.finditer(sources_match.group(1) if sources_match else "")
        }
        for line_no, line in enumerate(text.splitlines(), start=1):
            stripped = line.strip()
            spans = (
                [stripped]
                if re.match(r"^Evidence \(", stripped)
                else [
                    stripped[m.start() :]
                    for m in re.finditer(r"\(Ev:", _without_inline_code(stripped))
                ]
            )
            for span in spans:
                prev_code: Optional[str] = None
                for token, in_parens in _evidence_tokens(span):
                    head, sep, symbol = (part.strip() for part in token.partition("::"))
                    problem: Optional[str] = None
                    if any(
                        ch in head for ch in " *<{"
                    ) or not EVIDENCE_FILE_SUFFIX_RE.search(head):
                        if sep and in_parens and prev_code:
                            problem = _symbol_problem(token, prev_code)
                        elif sep:
                            problem = f"Evidence :: token must start with a source file path: {token}"
                    else:
                        matches = [
                            t for t in tracked if t == head or t.endswith("/" + head)
                        ]
                        matches = (
                            [t for t in matches if ROOT / t in sources]
                            or [t for t in matches if t == head]
                            or matches
                        )
                        prev_code = None
                        if len(matches) != 1:
                            problem = (
                                f"ambiguous Evidence file {head!r} matches {', '.join(matches[:3])}"
                                if matches
                                else f"Evidence file {head!r} is not a git-tracked file on disk"
                            )
                        else:
                            if matches[0].endswith(CODE_SUFFIXES):
                                prev_code = matches[0]
                            if sep:
                                problem = _symbol_problem(symbol, matches[0])
                    if problem:
                        _fail(failures, page, line_no, problem)


def _main() -> int:
    """
    Run the KB link, reachability, and citation checks.

    :return: Process exit status.
    """
    pages = sorted(p.resolve() for p in WIKI.rglob("*.md"))
    files = sorted({AGENTS.resolve(), SOURCES.resolve(), *pages})
    failures: List[str] = []
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
