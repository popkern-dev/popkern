"""CalculiX heat-transfer oracle (PH-065): solves a part's steady
conduction with quadratic tetrahedra (DC3D10) and writes the surface
nodes' heat flux magnitude by surface type, for comparison with
popkern-analysis's thermal solve.

    python3 heat_oracle.py <part.step> <setup.toml> <out dir> <mesh size> \
        <curvature points> <conductivity W/(m K)>

The setup is popkern-analysis's own thermal file: faces selected by
plane, `fixed_temperature` on faces and `heat_flux` given as a total
power spread uniformly over the faces; every other face is insulated.
Units are mm and W, so the conductivity goes to the deck in W/(mm K)
and CalculiX's flux (W/mm²) is written in W/m². The mesh and its
numbering checks are oracle.py's.
"""
import json
import math
import os
import pathlib
import subprocess
import sys
import tomllib

import gmsh

from oracle import C3D10_EDGES, GMSH_TO_C3D10, plane_faces


def solve(out, deck):
    """Runs CalculiX on `deck` in `out`."""
    (out / "job.inp").write_text(deck)
    env = dict(os.environ, OMP_NUM_THREADS=os.environ.get("OMP_NUM_THREADS", "8"))
    run = subprocess.run(["ccx", "-i", "job"], cwd=out, env=env, capture_output=True, text=True)
    (out / "ccx.log").write_text(run.stdout + run.stderr)
    if run.returncode != 0 or "*ERROR" in run.stdout:
        sys.exit("ccx failed; see ccx.log")


def main():
    step, setup_path, out, size, curvature, conductivity = sys.argv[1:7]
    out = pathlib.Path(out)
    out.mkdir(parents=True, exist_ok=True)
    setup = tomllib.loads(pathlib.Path(setup_path).read_text())
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.occ.importShapes(step)
    gmsh.model.occ.synchronize()
    if len(gmsh.model.getEntities(3)) != 1:
        sys.exit("one solid expected")
    bb = gmsh.model.getBoundingBox(-1, -1)
    tol = 1e-6 * max(bb[3] - bb[0], bb[4] - bb[1], bb[5] - bb[2])
    faces = {}
    for name, entry in setup["faces"].items():
        if "plane" not in entry:
            sys.exit(f"face {name}: only plane selectors are supported")
        faces[name] = plane_faces(entry["plane"], tol)
        if not faces[name]:
            sys.exit(f"face {name}: no surface matches")
    gmsh.option.setNumber("Mesh.ElementOrder", 2)
    gmsh.option.setNumber("Mesh.SecondOrderLinear", 0)
    gmsh.option.setNumber("Mesh.HighOrderOptimize", 2)
    gmsh.option.setNumber("Mesh.MeshSizeMax", float(size))
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", float(curvature))
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.model.mesh.generate(3)

    tags, coords, _ = gmsh.model.mesh.getNodes()
    xyz = {int(t): coords[3 * i: 3 * i + 3] for i, t in enumerate(tags)}
    types, _, conn = gmsh.model.mesh.getElements(3)
    if list(types) != [11]:
        sys.exit(f"unexpected volume element types {list(types)}")
    tets = [[int(v) for v in conn[0][10 * i: 10 * i + 10]] for i in range(len(conn[0]) // 10)]
    tets = [[t[k] for k in GMSH_TO_C3D10] for t in tets]
    worst = 0.0
    for t in tets:
        c = [xyz[v] for v in t]
        for k, (a, b) in enumerate(C3D10_EDGES):
            mid = [(c[a][i] + c[b][i]) / 2 for i in range(3)]
            worst = max(worst, math.dist(mid, c[4 + k]) / math.dist(c[a], c[b]))
    if worst > 0.25:
        sys.exit(f"edge-node numbering check failed ({worst})")

    def triangles(surface_tags):
        tris = []
        for s in surface_tags:
            ty, _, nodes = gmsh.model.mesh.getElements(2, s)
            for k, t in enumerate(ty):
                if t != 9:
                    sys.exit(f"unexpected surface element type {t}")
                n = nodes[k]
                tris += [[int(v) for v in n[6 * i: 6 * i + 6]] for i in range(len(n) // 6)]
        return tris

    def area(tri):
        a, b, c = (xyz[v] for v in tri[:3])
        u = [b[i] - a[i] for i in range(3)]
        w = [c[i] - a[i] for i in range(3)]
        n = [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]]
        return 0.5 * math.sqrt(sum(x * x for x in n))

    # Surface nodes by surface: its tag, type and bounding box.
    surfaces = []
    for _, s in gmsh.model.getEntities(2):
        nodes = {v for tri in triangles([s]) for v in tri}
        surfaces.append((s, gmsh.model.getType(2, s), gmsh.model.getBoundingBox(2, s), nodes))
    fixed = {}
    for entry in setup.get("fixed_temperature", []):
        for name in entry["faces"]:
            for tri in triangles(faces[name]):
                for v in tri:
                    fixed[v] = entry["temperature_c"]
    flux_in = {}
    for entry in setup.get("heat_flux", []):
        tris = triangles([t for name in entry["faces"] for t in faces[name]])
        density = entry["power_w"] / sum(area(t) for t in tris)
        # A 6-node triangle under uniform flux: nothing at the corners,
        # a third of its heat at each edge node.
        for tri in tris:
            for v in tri[3:]:
                flux_in[v] = flux_in.get(v, 0.0) + density * area(tri) / 3.0
    gmsh.finalize()

    deck = ["*NODE"]
    deck += [f"{t},{x[0]:.13e},{x[1]:.13e},{x[2]:.13e}" for t, x in sorted(xyz.items())]
    deck.append("*ELEMENT,TYPE=DC3D10,ELSET=EALL")
    deck += [f"{k + 1}," + ",".join(str(v) for v in t) for k, t in enumerate(tets)]
    deck += [
        "*MATERIAL,NAME=M", "*CONDUCTIVITY", f"{float(conductivity) * 1e-3:.13e}",
        "*SOLID SECTION,ELSET=EALL,MATERIAL=M",
        "*STEP", "*HEAT TRANSFER,STEADY STATE", "1.,1.", "*BOUNDARY",
    ]
    deck += [f"{v},11,11,{t:.13e}" for v, t in sorted(fixed.items())]
    deck.append("*CFLUX")
    deck += [f"{v},11,{q:.13e}" for v, q in sorted(flux_in.items())]
    deck += ["*NODE FILE", "NT", "*EL FILE", "HFL", "*END STEP"]
    deck = "\n".join(deck) + "\n"
    # A solved deck is reused: Gmsh meshes the same input identically.
    solved = (out / "job.inp").exists() and (out / "job.inp").read_text() == deck
    if not (solved and (out / "job.frd").exists()):
        solve(out, deck)
    temp, flux, block = {}, {}, None
    for line in (out / "job.frd").read_text().splitlines():
        if line.startswith(" -4"):
            block = line.split()[1]
        elif line.startswith(" -1") and block in ("NDTEMP", "FLUX"):
            node = int(line[3:13])
            vals = [float(line[13 + 12 * i: 25 + 12 * i]) for i in range((len(line) - 13) // 12)]
            (temp if block == "NDTEMP" else flux)[node] = vals
        elif line.startswith(" -3"):
            block = None
    if not temp or not flux:
        sys.exit("no NDTEMP or FLUX block in job.frd")

    def magnitude(v):
        return 1e6 * math.sqrt(sum(c * c for c in flux[v][:3]))

    with open(out / "surface.csv", "w") as f:
        f.write("x,y,z,surface,heat_flux\n")
        for tag, _, _, nodes in surfaces:
            for v in sorted(nodes):
                x = xyz[v]
                f.write(f"{x[0]:.9g},{x[1]:.9g},{x[2]:.9g},{tag},{magnitude(v):.9g}\n")
    peaks = []
    for tag, kind, box, nodes in surfaces:
        v = max(nodes, key=magnitude)
        peaks.append({"surface": tag, "type": kind, "box": list(box),
                      "heat_flux": magnitude(v), "at": list(xyz[v])})
    summary = {
        "nodes": len(xyz),
        "elements": len(tets),
        "fixed_nodes": len(fixed),
        "flux_in_w": sum(flux_in.values()),
        "max_temperature_c": max(t[0] for t in temp.values()),
        "peaks": peaks,
        "edge_node_check": worst,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
