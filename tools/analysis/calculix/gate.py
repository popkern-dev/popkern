"""PH-075: popkern-analysis against CalculiX as a gate.

    python3 gate.py <popkern-analysis binary> <work dir> [--record] [--tolerance-scale S]

For each case in gates.toml: CalculiX solves the part and setup
(oracle.py or heat_oracle.py; cached in the work dir by mesh), ours
solves it through the CLI, and the ratio ours / CalculiX of each
checked quantity must lie within `tolerance` of the recorded `ratio`.
Exits 1 when any check misses. `--record` prints the measured ratios
(to write into gates.toml) without judging; `--tolerance-scale` scales
every tolerance (0 shows the gate failing).

Structural checks: `displacement` (largest displacement) and `peak`
(von Mises peak inside the case's region of interest, clear of the
clamp's edge singularity). Thermal check: `bore_flux`, the heat-flux
peak on one CalculiX surface (by tag) against ours from `--refine`'s
submodel on that face.
"""
import csv
import json
import math
import pathlib
import subprocess
import sys
import tomllib

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parent.parent.parent
YOUNGS, POISSON = 3019.0, 0.36  # PLA (FDM), crates/popkern-analysis/src/material/presets.rs
CONDUCTIVITY = 167.0  # Aluminium 6061-T6


def inside(p, roi):
    return all(roi["lo"][i] <= p[i] <= roi["hi"][i] for i in range(3))


def ccx(case, work):
    """CalculiX's numbers for `case`, solving when not cached."""
    size, curvature = case["mesh"]
    out = work / case["name"] / f"ccx-{size}-{curvature}"
    step, setup = REPO / case["step"], REPO / case["setup"]
    if not (out / "summary.json").exists():
        if case["kind"] == "thermal":
            cmd = [sys.executable, str(HERE / "heat_oracle.py"), str(step), str(setup), str(out),
                   str(size), str(curvature), str(CONDUCTIVITY)]
        else:
            cmd = [sys.executable, str(HERE / "oracle.py"), str(step), str(setup), str(out),
                   str(size), str(curvature), str(YOUNGS), str(POISSON)]
        subprocess.run(cmd, check=True, capture_output=True)
    summary = json.loads((out / "summary.json").read_text())
    if case["kind"] == "thermal":
        peak = next(p for p in summary["peaks"] if p["surface"] == case["surface"])
        return {"bore_flux": peak["heat_flux"]}
    peak = 0.0
    with open(out / "surface.csv") as f:
        for r in csv.DictReader(f):
            p = [float(r["x"]), float(r["y"]), float(r["z"])]
            if inside(p, case["roi"]):
                peak = max(peak, float(r["von_mises"]))
    return {"displacement": summary["max_displacement"], "peak": peak}


def ours(binary, case, work):
    """Our numbers for `case` through the CLI."""
    out = work / case["name"] / "ours"
    out.mkdir(parents=True, exist_ok=True)
    text = (REPO / case["setup"]).read_text()
    if "cell_size" in case:
        text += f"\n[options]\ntarget_cell_size = {case['cell_size']}\n"
    setup = out / "setup.toml"
    setup.write_text(text)
    args = [binary, "run", str(REPO / case["step"]), str(setup), "--out", str(out)]
    args += ["--refine"] if case["kind"] == "thermal" else ["--samples"]
    subprocess.run(args, check=False, capture_output=True)
    report = json.loads((out / "setup.report.json").read_text())
    if report["status"] != "solved":
        raise SystemExit(f"{case['name']}: ours {report['status']}: {report.get('error')}")
    if case["kind"] == "thermal":
        subs = report["result"]["submodels"] or []
        if not subs:
            raise SystemExit(f"{case['name']}: no submodel refined the bore")
        return {"bore_flux": max(s["peak_heat_flux"] for s in subs if "peak_heat_flux" in s)}
    peak = 0.0
    with open(out / "setup.samples.csv") as f:
        for r in csv.DictReader(f):
            p = [float(r["x"]), float(r["y"]), float(r["z"])]
            if inside(p, case["roi"]):
                peak = max(peak, float(r["von_mises"]))
    return {"displacement": report["result"]["max_displacement_mm"], "peak": peak}


def main():
    args = sys.argv[1:]
    record = "--record" in args
    scale = 1.0
    if "--tolerance-scale" in args:
        scale = float(args[args.index("--tolerance-scale") + 1])
    binary, work = args[0], pathlib.Path(args[1])
    cases = tomllib.loads((HERE / "gates.toml").read_text())["case"]
    misses = 0
    for case in cases:
        theirs, mine = ccx(case, work), ours(binary, case, work)
        for name, check in case["checks"].items():
            ratio = mine[name] / theirs[name]
            if record:
                print(f"{case['name']:24} {name:12} ratio {ratio:.4f}")
                continue
            tol = scale * check["tolerance"]
            ok = abs(ratio - check["ratio"]) <= tol
            misses += not ok
            mark = "ok  " if ok else "MISS"
            print(f"{mark} {case['name']:24} {name:12} ratio {ratio:.4f}"
                  f" (recorded {check['ratio']:.4f} ± {tol:.4f})")
    if misses:
        print(f"{misses} check(s) missed")
        sys.exit(1)


if __name__ == "__main__":
    main()
