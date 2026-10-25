"""PH-105: every frozen gallery part through the CLI under six setups.

    python3 sweep.py <popkern-analysis binary> <work dir> [--cap SECONDS]

Static between faces 0 and 1 (plain, --refine-check, --refine), radiating
thermal with --refine, clamped modal and a radiating transient. Typed
declines are fine (a multi-solid part supported on one solid declines
FreeMotion). Exits 1 on an exit code other than 0 or 2, a missing
report, an Internal error kind, or a run longer than the cap (default
300 s; the sweep that made this tool found a 30-minute run, PH-104).
"""
import json
import pathlib
import subprocess
import sys
import time

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parent.parent.parent
RUNS = [
    ("static", []),
    ("static", ["--refine-check"]),
    ("static", ["--refine"]),
    ("thermal", ["--refine"]),
    ("modal", []),
    ("transient", []),
]


def main():
    args = sys.argv[1:]
    cap = 300.0
    if "--cap" in args:
        cap = float(args[args.index("--cap") + 1])
    binary, work = args[0], pathlib.Path(args[1])
    failures = []
    for step in sorted((REPO / "tests/fixtures/analysis/gallery").glob("*.step")):
        for setup, flags in RUNS:
            name = setup + "".join(flags).replace("--", "_")
            out = work / step.stem / name
            out.mkdir(parents=True, exist_ok=True)
            start = time.time()
            try:
                p = subprocess.run(
                    [binary, "run", str(step), str(HERE / "setups" / f"{setup}.toml"), "--out", str(out)] + flags,
                    capture_output=True, text=True, timeout=cap)
                code = p.returncode
            except subprocess.TimeoutExpired:
                code = None
            elapsed = time.time() - start
            report = out / f"{setup}.report.json"
            status, kind = "-", ""
            if code is not None and report.exists():
                r = json.loads(report.read_text())
                status = r.get("status", "-")
                kind = (r.get("error") or {}).get("kind", "")
            problem = None
            if code is None:
                problem = f"over the {cap:.0f} s cap"
            elif code not in (0, 2):
                problem = f"exit {code}: {p.stderr.strip()[:200]}"
            elif not report.exists():
                problem = "no report"
            elif kind == "Internal":
                problem = "internal error"
            mark = "FAIL" if problem else "ok  "
            print(f"{mark} {step.stem:22} {name:22} {status:9} {kind:20} {elapsed:6.1f}s"
                  + (f"  {problem}" if problem else ""), flush=True)
            if problem:
                failures.append(f"{step.stem} {name}: {problem}")
    if failures:
        print(f"{len(failures)} run(s) failed")
        sys.exit(1)


if __name__ == "__main__":
    main()
