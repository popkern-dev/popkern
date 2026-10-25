"""CalculiX oracle (PH-050): solves a part and setup with quadratic
tetrahedra and writes the surface nodes' von Mises and the largest
displacement, for comparison with popkern-analysis.

    python3 oracle.py <part.step> <setup.toml> <out dir> <mesh size> \
        <curvature points> <youngs MPa> <poisson>

Numbers go to the deck with 14 significant digits: CalculiX refuses
longer fields.

The setup is popkern-analysis's own file: faces selected by plane
(normal and a point), fixed supports, and loads given as a total force
spread uniformly over the faces. Gmsh meshes the STEP (OCC kernel,
second-order, curved), the script writes the CalculiX deck (C3D10,
consistent nodal loads of the 6-node face triangles) and reads the
.frd results (nodal stresses extrapolated and averaged by CalculiX).
"""
import json
import math
import os
import pathlib
import subprocess
import sys
import tomllib

import gmsh

# Gmsh's 10-node tetrahedron numbers its last two edge nodes (2-3, 1-3)
# the other way round from CalculiX's C3D10 (1-3, 2-3), 0-based.
GMSH_TO_C3D10 = [0, 1, 2, 3, 4, 5, 6, 7, 9, 8]
# Corner pairs of the C3D10 edge nodes 4..9.
C3D10_EDGES = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)]


def plane_faces(selector, tol):
    """Surface tags lying in the plane through `through` with outward
    normal `normal` (normal checked at the face's parametric centre)."""
    n = selector["normal"]
    length = math.sqrt(sum(c * c for c in n))
    n = [c / length for c in n]
    p0 = selector["through"]
    out = []
    for _, tag in gmsh.model.getEntities(2):
        if gmsh.model.getType(2, tag) != "Plane":
            continue
        lo, hi = gmsh.model.getParametrizationBounds(2, tag)
        uv = [0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1])]
        x = gmsh.model.getValue(2, tag, uv)
        normal = gmsh.model.getNormal(tag, uv)
        if sum(a * b for a, b in zip(normal, n)) < 1.0 - 1e-6:
            continue
        if abs(sum((x[i] - p0[i]) * n[i] for i in range(3))) <= tol:
            out.append(tag)
    return out


def main():
    step, setup_path, out, size, curvature, youngs, poisson = sys.argv[1:8]
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
    # Check the numbering: edge nodes near their corners' midpoints, and
    # positive corner volumes.
    worst = 0.0
    for t in tets:
        c = [xyz[v] for v in t]
        for k, (a, b) in enumerate(C3D10_EDGES):
            mid = [(c[a][i] + c[b][i]) / 2 for i in range(3)]
            edge = math.dist(c[a], c[b])
            worst = max(worst, math.dist(mid, c[4 + k]) / edge)
        u, v, w = ([c[j][i] - c[0][i] for i in range(3)] for j in (1, 2, 3))
        vol = (u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0])
               + u[2] * (v[0] * w[1] - v[1] * w[0]))
        if vol <= 0:
            sys.exit("non-positive tetrahedron")
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

    surface = sorted({v for tri in triangles([t for _, t in gmsh.model.getEntities(2)]) for v in tri})
    fixed = set()
    for support in setup.get("support", []):
        if support.get("kind") != "fixed":
            sys.exit("only fixed supports are supported")
        for name in support["faces"]:
            fixed |= {v for tri in triangles(faces[name]) for v in tri}
    loads = {}
    for load in setup.get("load", []):
        tris = triangles([t for name in load["faces"] for t in faces[name]])
        area = lambda tri: 0.5 * math.sqrt(sum(
            x * x for x in [
                (xyz[tri[1]][1] - xyz[tri[0]][1]) * (xyz[tri[2]][2] - xyz[tri[0]][2])
                - (xyz[tri[1]][2] - xyz[tri[0]][2]) * (xyz[tri[2]][1] - xyz[tri[0]][1]),
                (xyz[tri[1]][2] - xyz[tri[0]][2]) * (xyz[tri[2]][0] - xyz[tri[0]][0])
                - (xyz[tri[1]][0] - xyz[tri[0]][0]) * (xyz[tri[2]][2] - xyz[tri[0]][2]),
                (xyz[tri[1]][0] - xyz[tri[0]][0]) * (xyz[tri[2]][1] - xyz[tri[0]][1])
                - (xyz[tri[1]][1] - xyz[tri[0]][1]) * (xyz[tri[2]][0] - xyz[tri[0]][0]),
            ]))
        total = sum(area(t) for t in tris)
        d = load["direction"]
        length = math.sqrt(sum(c * c for c in d))
        traction = [load["force_n"] * c / length / total for c in d]
        # A 6-node triangle under uniform traction: nothing at the
        # corners, a third of its force at each edge node.
        for tri in tris:
            a = area(tri)
            for v in tri[3:]:
                f = loads.setdefault(v, [0.0, 0.0, 0.0])
                for i in range(3):
                    f[i] += traction[i] * a / 3.0
    gmsh.finalize()

    deck = ["*NODE"]
    deck += [f"{t},{x[0]:.13e},{x[1]:.13e},{x[2]:.13e}" for t, x in sorted(xyz.items())]
    deck.append("*ELEMENT,TYPE=C3D10,ELSET=EALL")
    deck += [f"{k + 1}," + ",".join(str(v) for v in t) for k, t in enumerate(tets)]
    deck.append("*NSET,NSET=FIXED")
    deck += [str(v) for v in sorted(fixed)]
    deck += [
        "*MATERIAL,NAME=M", "*ELASTIC", f"{youngs},{poisson}",
        "*SOLID SECTION,ELSET=EALL,MATERIAL=M",
        "*STEP", "*STATIC", "*BOUNDARY", "FIXED,1,3", "*CLOAD",
    ]
    for v, f in sorted(loads.items()):
        deck += [f"{v},{i + 1},{f[i]:.13e}" for i in range(3) if f[i] != 0.0]
    deck += ["*NODE FILE", "U", "*EL FILE", "S", "*END STEP"]
    (out / "job.inp").write_text("\n".join(deck) + "\n")
    env = dict(os.environ, OMP_NUM_THREADS=os.environ.get("OMP_NUM_THREADS", "8"))
    run = subprocess.run(["ccx", "-i", "job"], cwd=out, env=env, capture_output=True, text=True)
    (out / "ccx.log").write_text(run.stdout + run.stderr)
    if run.returncode != 0 or "*ERROR" in run.stdout:
        sys.exit("ccx failed; see ccx.log")

    disp, stress, block = {}, {}, None
    for line in (out / "job.frd").read_text().splitlines():
        if line.startswith(" -4"):
            block = line.split()[1]
        elif line.startswith(" -1") and block in ("DISP", "STRESS"):
            node = int(line[3:13])
            vals = [float(line[13 + 12 * i: 25 + 12 * i]) for i in range((len(line) - 13) // 12)]
            (disp if block == "DISP" else stress)[node] = vals
        elif line.startswith(" -3"):
            block = None

    def von_mises(s):
        xx, yy, zz, xy, yz, zx = s[:6]
        return math.sqrt(0.5 * ((xx - yy) ** 2 + (yy - zz) ** 2 + (zz - xx) ** 2)
                         + 3.0 * (xy * xy + yz * yz + zx * zx))

    max_disp = max(math.sqrt(sum(c * c for c in u[:3])) for u in disp.values())
    with open(out / "surface.csv", "w") as f:
        f.write("x,y,z,von_mises\n")
        for v in surface:
            x = xyz[v]
            f.write(f"{x[0]:.9g},{x[1]:.9g},{x[2]:.9g},{von_mises(stress[v]):.9g}\n")
    summary = {
        "nodes": len(xyz),
        "elements": len(tets),
        "surface_nodes": len(surface),
        "fixed_nodes": len(fixed),
        "applied": [sum(f[i] for f in loads.values()) for i in range(3)],
        "max_displacement": max_disp,
        "edge_node_check": worst,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
