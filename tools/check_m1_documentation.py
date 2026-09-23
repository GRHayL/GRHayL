#!/usr/bin/env python3
"""Check repository-relative links in the maintained Radiation M1 documents."""

from pathlib import Path
import re
import sys


ROOT = Path(__file__).resolve().parents[1]
PAGE_ROOTS = (
    ROOT / "GRHayL/Radiation/M1_INTEGRATION_CONTRACT.md",
    ROOT / "Unit_Tests/README.m1.md",
    ROOT / "Unit_Tests/data/m1_thcm1",
    ROOT / "wiki/gems/radiation-m1",
)
LINK = re.compile(r"\[[^\]]*\]\(([^)\n]+)\)")
FENCE = re.compile(r"^\s*(`{3,}|~{3,})")
SCHEME = re.compile(r"^[A-Za-z][A-Za-z0-9+.-]*:")


def pages() -> list[Path]:
    result = []
    for entry in PAGE_ROOTS:
        if entry.is_file():
            result.append(entry)
        elif entry.is_dir():
            result.extend(sorted(entry.rglob("*.md")))
        else:
            print(f"missing documentation root: {entry}", file=sys.stderr)
    return result


def local_target(raw_target: str) -> str | None:
    target = raw_target.strip()
    if target.startswith("<"):
        end = target.find(">")
        if end < 0:
            return None
        target = target[1:end]
    else:
        target = target.split(None, 1)[0]
    target = target.split("#", 1)[0]
    if not target or target.startswith("//") or SCHEME.match(target):
        return None
    return target


def main() -> int:
    failures = 0
    for page in pages():
        in_fence = False
        fence_marker = None
        for line_number, line in enumerate(page.read_text(encoding="utf-8").splitlines(), 1):
            fence = FENCE.match(line)
            if fence:
                marker = fence.group(1)[0]
                if not in_fence:
                    in_fence = True
                    fence_marker = marker
                elif marker == fence_marker:
                    in_fence = False
                    fence_marker = None
                continue
            if in_fence:
                continue
            for raw_target in LINK.findall(line):
                target = local_target(raw_target)
                if target is None:
                    continue
                candidate = (page.parent / target).resolve()
                try:
                    candidate.relative_to(ROOT)
                except ValueError:
                    print(f"{page.relative_to(ROOT)}:{line_number}: escaping target: {raw_target}")
                    failures += 1
                    continue
                if not candidate.exists():
                    print(f"{page.relative_to(ROOT)}:{line_number}: missing target: {raw_target}")
                    failures += 1
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
