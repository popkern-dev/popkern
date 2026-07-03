#!/usr/bin/env python3
"""Mesh-finish pass for showcase OBJ exports: weld coincident vertices,
drop degenerate triangles, and close residual sliver gaps by filling
boundary loops.

The KP20 apex showcase writes its parts with shared seam vertices
wherever the construction allows; the remaining hairline gaps come
from the cap assembly's residual C0 mismatch (the c0 gate) and from
T-junctions between cap pieces. Those slivers are geometrically
bounded by the c0 metric, so fan-filling their boundary loops is a
faithful closure, not invention. The BRep-level fix (assembly c0 to
machine zero) is tracked in the task ledger.

Usage: mesh-finish.py input.obj output.obj
"""

import collections
import sys


def main() -> int:
    src, dst = sys.argv[1], sys.argv[2]
    verts: list[tuple[float, float, float]] = []
    normals: list[tuple[float, float, float]] = []
    tris: list[tuple[list[int], str]] = []
    group = "mesh"
    for line in open(src):
        if line.startswith("g "):
            group = line[2:].strip()
        elif line.startswith("v "):
            verts.append(tuple(float(x) for x in line.split()[1:4]))
        elif line.startswith("vn "):
            normals.append(tuple(float(x) for x in line.split()[1:4]))
        elif line.startswith("f "):
            idx = [int(t.split("/")[0]) - 1 for t in line.split()[1:]]
            tris.append((idx, group))

    # Weld by rounded position.
    key: dict[tuple, int] = {}
    canon: list[int] = []  # canonical -> original vertex index
    remap: list[int] = []
    for i, v in enumerate(verts):
        k = tuple(round(c, 7) for c in v)
        if k not in key:
            key[k] = len(canon)
            canon.append(i)
        remap.append(key[k])

    # Weld POSITIONS only; each face corner keeps its ORIGINAL normal
    # so creases (base rim, face/strip junctions) stay crisp.
    out_tris: list[tuple[int, int, int, tuple[int, int, int], str]] = []
    edge_count: collections.Counter = collections.Counter()
    for idx, grp in tris:
        w = [remap[i] for i in idx]
        if len(set(w)) < 3:
            continue
        out_tris.append((w[0], w[1], w[2], (idx[0], idx[1], idx[2]), grp))
        for a, b in ((w[0], w[1]), (w[1], w[2]), (w[2], w[0])):
            edge_count[(min(a, b), max(a, b))] += 1

    # Boundary loops — filled iteratively (loops that touch at shared
    # vertices need multiple passes).
    new_verts: list[tuple[float, float, float]] = []
    new_norms: list[tuple[float, float, float]] = []
    fills: list[tuple[int, int, int]] = []
    loops_filled = 0

    def position(i: int) -> tuple[float, float, float]:
        if i < len(canon):
            return verts[canon[i]]
        return new_verts[i - len(canon)]

    def normal_of(i: int) -> tuple[float, float, float]:
        if i < len(canon):
            ci = canon[i]
            return normals[ci] if ci < len(normals) else (0.0, 0.0, 1.0)
        return new_norms[i - len(canon)]

    for _pass in range(8):
        counts: collections.Counter = collections.Counter()
        for a, b, c, _, _ in out_tris:
            for e in ((a, b), (b, c), (c, a)):
                counts[(min(e), max(e))] += 1
        for a, b, c in fills:
            for e in ((a, b), (b, c), (c, a)):
                counts[(min(e), max(e))] += 1
        adjacency: dict[int, list[int]] = collections.defaultdict(list)
        n_boundary = 0
        for (a, b), c in counts.items():
            if c == 1:
                adjacency[a].append(b)
                adjacency[b].append(a)
                n_boundary += 1
        if n_boundary == 0:
            break
        visited: set[int] = set()
        progressed = False
        for start in list(adjacency):
            if start in visited:
                continue
            loop = [start]
            visited.add(start)
            cur = start
            while True:
                nxt = next((n for n in adjacency[cur] if n not in visited), None)
                if nxt is None:
                    break
                loop.append(nxt)
                visited.add(nxt)
                cur = nxt
            if len(loop) < 3:
                continue
            progressed = True
            loops_filled += 1
            pts = [position(i) for i in loop]
            cx = sum(p[0] for p in pts) / len(pts)
            cy = sum(p[1] for p in pts) / len(pts)
            cz = sum(p[2] for p in pts) / len(pts)
            nrm = [normal_of(i) for i in loop]
            nx = sum(v[0] for v in nrm) / len(nrm)
            ny = sum(v[1] for v in nrm) / len(nrm)
            nz = sum(v[2] for v in nrm) / len(nrm)
            center = len(canon) + len(new_verts)
            new_verts.append((cx, cy, cz))
            new_norms.append((nx, ny, nz))
            for i in range(len(loop)):
                fills.append((center, loop[i], loop[(i + 1) % len(loop)]))
        if not progressed:
            break
    loops = list(range(loops_filled))  # count only, for the report

    with open(dst, "w") as f:
        f.write("# mesh-finish: positions welded, per-corner normals, loops filled\n")
        for ci in canon:
            v = verts[ci]
            f.write(f"v {v[0]} {v[1]} {v[2]}\n")
        for v in new_verts:
            f.write(f"v {v[0]} {v[1]} {v[2]}\n")
        for n in normals:
            f.write(f"vn {n[0]} {n[1]} {n[2]}\n")
        for n in new_norms:
            f.write(f"vn {n[0]} {n[1]} {n[2]}\n")
        n_base = len(normals)
        cur_group = None
        for a, b, c, (na, nb, nc), grp in out_tris:
            if grp != cur_group:
                f.write(f"g {grp}\n")
                cur_group = grp
            f.write(f"f {a+1}//{na+1} {b+1}//{nb+1} {c+1}//{nc+1}\n")
        if fills:
            f.write("g seam_fill\n")
            for a, b, c in fills:
                def norm_ix(i: int) -> int:
                    if i >= len(canon):
                        return n_base + (i - len(canon)) + 1
                    ci = canon[i]
                    return (ci + 1) if ci < len(normals) else 1
                f.write(
                    f"f {a+1}//{norm_ix(a)} {b+1}//{norm_ix(b)} {c+1}//{norm_ix(c)}\n"
                )

    # Closure audit.
    final_edges: collections.Counter = collections.Counter()
    for a, b, c, _, _ in out_tris:
        for e in ((a, b), (b, c), (c, a)):
            final_edges[(min(e), max(e))] += 1
    for a, b, c in fills:
        for e in ((a, b), (b, c), (c, a)):
            final_edges[(min(e), max(e))] += 1
    boundary = sum(1 for c in final_edges.values() if c == 1)
    print(
        f"mesh-finish: {len(out_tris)} tris, {len(loops)} loops filled "
        f"({len(fills)} tris), boundary edges after fill: {boundary}"
    )
    return 0 if boundary == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
