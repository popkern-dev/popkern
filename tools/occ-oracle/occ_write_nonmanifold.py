#!/usr/bin/env python3
"""Write non-manifold STEP fixtures with OpenCASCADE (PK-980).

Each fixture is one MANIFOLD_SOLID_BREP whose CLOSED_SHELL holds the
faces of several solids that touch: the solids are made to share their
contact edges and vertices by a general fuse, and their faces go into
one shell, so the file lists a contact edge with four face uses (or a
contact vertex shared by two closed sets of faces). OCC reads such a
file back as one solid; the expected pieces are the input solids.

Usage:
    <ocp-venv>/bin/python tools/occ-oracle/occ_write_nonmanifold.py OUT_DIR

Writes OUT_DIR/<name>.step for each fixture and prints one JSON object:
    {name: {"pieces": [volume, ...ascending], "occ_solids": [volume, ...]}}
where occ_solids are the solids OCC reads back from the written file.
"""

import json
import os
import sys


def box(x0, y0, z0, x1, y1, z1):
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
    from OCP.gp import gp_Pnt

    return BRepPrimAPI_MakeBox(gp_Pnt(x0, y0, z0), gp_Pnt(x1, y1, z1)).Shape()


def cylinder(x, y, z, r, h, axis=(0, 0, 1)):
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeCylinder
    from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt

    return BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(x, y, z), gp_Dir(*axis)), r, h).Shape()


def fuse_all(*shapes):
    """One solid: the shapes fused (a solid that may touch itself)."""
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse

    out = shapes[0]
    for s in shapes[1:]:
        out = BRepAlgoAPI_Fuse(out, s).Shape()
    return out


FIXTURES = {
    # Two cubes sharing one edge (x = y = 1, along z).
    "nm_cubes_on_an_edge": lambda: [box(0, 0, 0, 1, 1, 1), box(1, 1, 0, 2, 2, 1)],
    # Two cubes sharing one corner vertex.
    "nm_cubes_at_a_corner": lambda: [box(0, 0, 0, 1, 1, 1), box(1, 1, 1, 2, 2, 2)],
    # Three blocks round one edge, at mm scale (x = y = 10, along z).
    "nm_three_blocks_on_an_edge": lambda: [
        box(0, 0, 0, 10, 10, 5),
        box(10, 10, 0, 20, 20, 5),
        box(0, 10, -5, 10, 20, 0),
    ],
    # A cylinder lying on a block's top face along a line.
    "nm_cylinder_on_a_block": lambda: [
        box(-3, -3, 0, 3, 3, 2),
        cylinder(-3, 0, 3, 1, 6, axis=(1, 0, 0)),
    ],
    # One solid touching itself along an edge (PK-1070): two cubes on an
    # edge joined by a slab over both, fused, so the cubes' contact edge
    # has four face uses within one closed piece.
    "nm_one_piece_touching_itself": lambda: [
        fuse_all(
            box(0, 0, 0, 1, 1, 1),
            box(1, 1, 0, 2, 2, 1),
            box(0, 0, 1, 2, 2, 1.5),
        )
    ],
    # Two boxes sharing a face, the face listed once (PK-1071): an
    # internal face whose edges have three face uses.
    "nm_boxes_sharing_a_face": lambda: [box(0, 0, 0, 1, 1, 1), box(0, 0, 1, 1, 1, 2)],
}


def cut_all(shape, *tools):
    """`shape` with every tool cut away."""
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut

    for t in tools:
        shape = BRepAlgoAPI_Cut(shape, t).Shape()
    return shape


# Solids written as OpenCASCADE builds them (no one-shell assembly).
DIRECT_FIXTURES = {
    # A block with two cavities touching along an edge (PK-1072): a
    # BREP_WITH_VOIDS whose void shells share that edge.
    "nm_touching_cavities": lambda: cut_all(
        box(0, 0, 0, 4, 4, 4), box(1, 1, 1, 2, 2, 2), box(2, 2, 1, 3, 3, 2)
    ),
}


def fused_pieces(solids):
    """The input solids sharing their contacts: a general fuse result."""
    from OCP.BOPAlgo import BOPAlgo_Builder
    from OCP.TopAbs import TopAbs_SOLID
    from OCP.TopExp import TopExp_Explorer

    builder = BOPAlgo_Builder()
    for s in solids:
        builder.AddArgument(s)
    builder.Perform()
    if builder.HasErrors():
        raise SystemExit("general fuse failed")
    out = []
    explorer = TopExp_Explorer(builder.Shape(), TopAbs_SOLID)
    while explorer.More():
        out.append(explorer.Current())
        explorer.Next()
    return out


def one_shell_solid(pieces):
    """One solid whose single shell holds every piece's faces."""
    from OCP.BRep import BRep_Builder
    from OCP.TopAbs import TopAbs_FACE
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS, TopoDS_Shell, TopoDS_Solid

    builder = BRep_Builder()
    shell = TopoDS_Shell()
    builder.MakeShell(shell)
    added = []
    for piece in pieces:
        faces = TopExp_Explorer(piece, TopAbs_FACE)
        while faces.More():
            face = TopoDS.Face_s(faces.Current())
            # A face two pieces share (an internal face) is listed once.
            if not any(face.IsSame(f) for f in added):
                builder.Add(shell, face)
                added.append(face)
            faces.Next()
    shell.Closed(True)
    solid = TopoDS_Solid()
    builder.MakeSolid(solid)
    builder.Add(solid, shell)
    return solid


def volume(shape):
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps

    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, props, 1e-9)
    return props.Mass()


def read_back(path):
    from OCP.STEPControl import STEPControl_Reader
    from OCP.TopAbs import TopAbs_SOLID
    from OCP.TopExp import TopExp_Explorer

    reader = STEPControl_Reader()
    if reader.ReadFile(path) != 1:
        return None
    reader.TransferRoots()
    out = []
    for i in range(1, reader.NbShapes() + 1):
        explorer = TopExp_Explorer(reader.Shape(i), TopAbs_SOLID)
        while explorer.More():
            out.append(volume(explorer.Current()))
            explorer.Next()
    return sorted(out)


def write(shape, path):
    from OCP.IFSelect import IFSelect_RetDone
    from OCP.STEPControl import STEPControl_AsIs, STEPControl_Writer

    writer = STEPControl_Writer()
    if writer.Transfer(shape, STEPControl_AsIs) != IFSelect_RetDone:
        raise SystemExit(f"transfer failed for {path}")
    if writer.Write(path) != IFSelect_RetDone:
        raise SystemExit(f"write failed for {path}")


def main(out_dir):
    os.makedirs(out_dir, exist_ok=True)
    report = {}
    for name, make in FIXTURES.items():
        solids = make()
        # One solid (already fused) needs no general fuse.
        pieces = fused_pieces(solids) if len(solids) > 1 else solids
        path = os.path.join(out_dir, f"{name}.step")
        write(one_shell_solid(pieces), path)
        report[name] = {
            "pieces": sorted(volume(p) for p in pieces),
            "occ_solids": read_back(path),
        }
    for name, make in DIRECT_FIXTURES.items():
        shape = make()
        path = os.path.join(out_dir, f"{name}.step")
        write(shape, path)
        report[name] = {"pieces": [volume(shape)], "occ_solids": read_back(path)}
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main(sys.argv[1])
