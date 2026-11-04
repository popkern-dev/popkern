#!/usr/bin/env python3
"""OpenCASCADE boolean oracle: an independent union (or cut) of popkern's
operands (PK-758, PK-751).

Usage:
    <ocp-venv>/bin/python tools/occ-oracle/occ_fuse.py [--op fuse|cut|common|exactly-one] A.step B.step [C.step ...]

Reads every solid of every file and folds them left to right with
BRepAlgoAPI_Fuse (default), BRepAlgoAPI_Cut (the first operand minus
each later one) or BRepAlgoAPI_Common. exactly-one (the N-body
symmetric difference, PK-1048) is the compound of each operand cut by
the fuse of the others: the pieces are disjoint, so its volume is
their sum (OCC has no N-body exactly-one of its own). It then reports the result's validity and adaptive volume
(1e-9 relative, BRepGProp). A popkern result served from the same
operands is compared against this volume: the two kernels share no
boolean code.

Output: one JSON object on stdout:
    done      the fuse completed
    valid     BRepCheck_Analyzer on the fused shape
    solids    solid count of the fused shape
    volume    adaptive volume (mm^3)
"""

import json
import sys


def load(path):
    from OCP.STEPControl import STEPControl_Reader

    reader = STEPControl_Reader()
    if reader.ReadFile(path) != 1:
        raise SystemExit(json.dumps({"done": False, "error": f"cannot read {path}"}))
    reader.TransferRoots()
    return reader.OneShape()


def main(paths, op="fuse"):
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
    from OCP.BRepCheck import BRepCheck_Analyzer
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    from OCP.TopAbs import TopAbs_SOLID
    from OCP.TopExp import TopExp_Explorer

    from OCP.BRep import BRep_Builder
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Common
    from OCP.TopoDS import TopoDS_Compound

    def fold(algo_of, first, rest):
        acc = first
        for other in rest:
            algo = algo_of(acc, other)
            algo.Build()
            if not algo.IsDone():
                return None
            acc = algo.Shape()
        return acc

    shapes = [load(p) for p in paths]
    algos = {"cut": BRepAlgoAPI_Cut, "common": BRepAlgoAPI_Common, "fuse": BRepAlgoAPI_Fuse}
    if op == "exactly-one":
        compound = TopoDS_Compound()
        builder = BRep_Builder()
        builder.MakeCompound(compound)
        fused = compound
        for i, shape in enumerate(shapes):
            others = shapes[:i] + shapes[i + 1 :]
            rest = fold(BRepAlgoAPI_Fuse, others[0], others[1:])
            piece = None if rest is None else fold(BRepAlgoAPI_Cut, shape, [rest])
            if piece is None:
                print(json.dumps({"done": False}))
                return 1
            builder.Add(compound, piece)
    else:
        fused = fold(algos[op], shapes[0], shapes[1:])
        if fused is None:
            print(json.dumps({"done": False}))
            return 1
    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(fused, props, 1e-9)
    solids = 0
    explorer = TopExp_Explorer(fused, TopAbs_SOLID)
    while explorer.More():
        solids += 1
        explorer.Next()
    print(
        json.dumps(
            {
                "done": True,
                "valid": BRepCheck_Analyzer(fused).IsValid(),
                "solids": solids,
                "volume": abs(props.Mass()),
            }
        )
    )
    return 0


if __name__ == "__main__":
    args = sys.argv[1:]
    op = "fuse"
    if len(args) >= 2 and args[0] == "--op":
        op, args = args[1], args[2:]
    sys.exit(main(args, op))
