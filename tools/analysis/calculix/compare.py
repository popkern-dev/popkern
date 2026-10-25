"""PH-050: popkern-analysis against the CalculiX oracle.

    python3 compare.py <popkern-analysis binary> <work dir> [case ...]

For each case in cases.toml: our solve at the default cell budget and
at each listed cell size (surface samples via --samples), CalculiX at
each listed mesh (oracle.py; reused when its summary exists), then the
ratios ours / CalculiX's finest mesh of the largest displacement and of
the peak von Mises inside the case's region of interest, which keeps
clear of the clamp's edge singularity. Prints a table and writes
results.json in the work dir.
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


def inside(p, roi):
    return all(roi["lo"][i] <= p[i] <= roi["hi"][i] for i in range(3))


def roi_peak(rows, roi):
    best = (0.0, None)
    for p, vm in rows:
        if inside(p, roi) and vm > best[0]:
            best = (vm, p)
    return best


def ours(binary, case, work, cell_size):
    tag = "default" if cell_size is None else f"h{cell_size}"
    out = work / case["part"] / f"ours-{tag}"
    out.mkdir(parents=True, exist_ok=True)
    setup = (HERE / "cases" / f"{case['part']}.setup.toml").read_text()
    if cell_size is not None:
        setup += f"\n[options]\ntarget_cell_size = {cell_size}\nmax_active_cells = 2000000\n"
    path = out / "setup.toml"
    path.write_text(setup)
    step = REPO / "tests" / "fixtures" / "analysis" / "gallery" / f"{case['part']}.step"
    subprocess.run([binary, "run", str(step), str(path), "--out", str(out), "--samples"],
                   check=False, capture_output=True)
    report = json.loads((out / "setup.report.json").read_text())
    if report["status"] != "solved":
        return {"status": report["status"], "error": report.get("error")}
    with open(out / "setup.samples.csv") as f:
        rows = [([float(r["x"]), float(r["y"]), float(r["z"])], float(r["von_mises"]))
                for r in csv.DictReader(f)]
    h = report["solve"]["cell_size"]
    peak, at = roi_peak(rows, case["roi"])
    return {
        "status": "solved",
        "cell_size": h,
        "cells_active": report["solve"]["cells_active"],
        "feature_cells": case["feature_radius"] / h,
        "max_displacement": report["result"]["max_displacement_mm"],
        "roi_peak": peak,
        "roi_peak_at": at,
    }


def oracle(case, work, size, curvature):
    out = work / case["part"] / f"ccx-{size}-{curvature}"
    if not (out / "summary.json").exists():
        step = REPO / "tests" / "fixtures" / "analysis" / "gallery" / f"{case['part']}.step"
        setup = HERE / "cases" / f"{case['part']}.setup.toml"
        subprocess.run([sys.executable, str(HERE / "oracle.py"), str(step), str(setup), str(out),
                        str(size), str(curvature), str(YOUNGS), str(POISSON)], check=True)
    summary = json.loads((out / "summary.json").read_text())
    with open(out / "surface.csv") as f:
        rows = [([float(r["x"]), float(r["y"]), float(r["z"])], float(r["von_mises"]))
                for r in csv.DictReader(f)]
    peak, at = roi_peak(rows, case["roi"])
    return {
        "mesh": [size, curvature],
        "nodes": summary["nodes"],
        "max_displacement": summary["max_displacement"],
        "roi_peak": peak,
        "roi_peak_at": at,
    }


def main():
    binary, work = sys.argv[1], pathlib.Path(sys.argv[2])
    wanted = set(sys.argv[3:])
    cases = tomllib.loads((HERE / "cases.toml").read_text())["case"]
    results = []
    for case in cases:
        if wanted and case["part"] not in wanted:
            continue
        ccx = [oracle(case, work, s, c) for s, c in case["meshes"]]
        ref = ccx[-1]
        runs = [ours(binary, case, work, None)]
        runs += [ours(binary, case, work, h) for h in case.get("cell_sizes", [])]
        for r in runs:
            if r["status"] == "solved":
                r["displacement_ratio"] = r["max_displacement"] / ref["max_displacement"]
                r["peak_ratio"] = r["roi_peak"] / ref["roi_peak"]
                r["peak_offset"] = math.dist(r["roi_peak_at"], ref["roi_peak_at"])
        results.append({"part": case["part"], "load": case["load"], "ccx": ccx, "ours": runs})
        print(f"\n{case['part']} ({case['load']})")
        for c in ccx:
            print(f"  ccx mesh {c['mesh']}: {c['nodes']} nodes, max |u| {c['max_displacement']:.6g},"
                  f" roi peak {c['roi_peak']:.6g} at {[round(x, 3) for x in c['roi_peak_at']]}")
        for r in runs:
            if r["status"] != "solved":
                print(f"  ours: {r}")
                continue
            print(f"  ours h {r['cell_size']:.4g} ({r['cells_active']} cells,"
                  f" {r['feature_cells']:.1f} cells/feature radius): u ratio"
                  f" {r['displacement_ratio']:.4f}, peak ratio {r['peak_ratio']:.4f}"
                  f" (offset {r['peak_offset']:.3g})")
    (work / "results.json").write_text(json.dumps(results, indent=1) + "\n")


if __name__ == "__main__":
    main()
