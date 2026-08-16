#!/usr/bin/env python3
"""Ranking-identity comparison for two kp20_bi6 battery audit runs.

Compares per-test captured output between two JUnit XML files produced
by `cargo nextest run -p popkern-core --profile battery-audit`
(.config/nextest.toml). The battery's ranking sweeps report candidate
orderings, verdict labels, and counts via eprintln; a solver change
that is mathematically identical but floating-point-different (PK-291
faer LU swap) may legitimately move printed floats by ULPs while every
ranking and verdict must stay identical.

Normalization: finite float literals are replaced by a placeholder
token (their exact digits are allowed to drift); `nan`/`inf` become
distinct tokens because finiteness is verdict-like; integers, names,
orderings, and all other text must match exactly.

Usage: battery-ranking-diff.py BASELINE.xml CANDIDATE.xml
Exit 0 when every test's normalized output matches and the test sets
are identical; exit 1 otherwise, with a unified diff per mismatch.

Reported mismatches still need triage, not blind acceptance or
rejection: argmax/argmin labels among tied candidates (residuals
identical to ~11 digits), hover counts at internal gates (candidates
measured ON a threshold), and κ-amplified measures on rejected-class
candidates flip legitimately under any summation reorder — the
PK-291 adoption audit measured all three classes. A changed label
between candidates whose printed measures genuinely separate, or any
change echoed by a test assertion failure, is a real regression.
"""

import difflib
import re
import sys
import xml.etree.ElementTree as ET

FLOAT = re.compile(
    r"(?<![\w.])[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)(?![\w.])"
    r"|(?<![\w.])[-+]?(?:\d+\.\d*|\.\d+)(?![\w.])"
)
NAN = re.compile(r"(?<![\w.])[-+]?(?:NaN|nan)(?![\w.])")
INF = re.compile(r"(?<![\w.])[-+]?(?:inf|Inf|INF|infinity|Infinity)(?![\w.])")


HARNESS_NOISE = re.compile(
    r"^(running \d+ tests?$|test result:|test \S+ \.\.\. "
    r"|test \S+ has been running for over )"
)


def normalize(text: str) -> str:
    # Drop libtest harness boilerplate: its filtered/total counts
    # change whenever unrelated tests are added to the binary, which
    # is not a ranking difference.
    lines = [
        line
        for line in text.splitlines()
        if line.strip() and not HARNESS_NOISE.match(line)
    ]
    text = "\n".join(lines)
    text = NAN.sub("<nan>", text)
    text = INF.sub("<inf>", text)
    return FLOAT.sub("<f>", text)


def load(path: str) -> dict[str, str]:
    tests = {}
    root = ET.parse(path).getroot()
    for case in root.iter("testcase"):
        key = f"{case.get('classname', '')}::{case.get('name', '')}"
        chunks = []
        for tag in ("system-out", "system-err"):
            node = case.find(tag)
            if node is not None and node.text:
                chunks.append(node.text)
        tests[key] = normalize("\n".join(chunks))
    return tests


def main() -> int:
    if len(sys.argv) != 3:
        print(__doc__.strip(), file=sys.stderr)
        return 2
    baseline = load(sys.argv[1])
    candidate = load(sys.argv[2])

    failures = 0
    only_base = sorted(set(baseline) - set(candidate))
    only_cand = sorted(set(candidate) - set(baseline))
    for name in only_base:
        print(f"MISSING in candidate: {name}")
        failures += 1
    for name in only_cand:
        print(f"EXTRA in candidate: {name}")
        failures += 1

    mismatched = 0
    for name in sorted(set(baseline) & set(candidate)):
        if baseline[name] != candidate[name]:
            mismatched += 1
            failures += 1
            print(f"\nOUTPUT MISMATCH: {name}")
            diff = difflib.unified_diff(
                baseline[name].splitlines(),
                candidate[name].splitlines(),
                fromfile="baseline",
                tofile="candidate",
                lineterm="",
                n=1,
            )
            for line in list(diff)[:60]:
                print(f"  {line}")

    total = len(set(baseline) | set(candidate))
    print(
        f"\n{total} tests compared: {mismatched} output mismatches, "
        f"{len(only_base)} missing, {len(only_cand)} extra"
    )
    if failures == 0:
        print("RANKING IDENTITY: PASS")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
