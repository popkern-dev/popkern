#!/usr/bin/env python3
"""OpenCASCADE fuse oracle: an independent union of popkern's operands
(PK-758).

Usage:
    <ocp-venv>/bin/python tools/occ-oracle/occ_fuse.py A.step B.step [C.step ...]

Reads every solid of every file, fuses them all with BRepAlgoAPI_Fuse,
and reports the fused shape's validity and adaptive volume (1e-9
relative, BRepGProp). A popkern union served from the same operands is
compared against this volume: the two kernels share no boolean code.

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


def main(paths):
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse
    from OCP.BRepCheck import BRepCheck_Analyzer
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    from OCP.TopAbs import TopAbs_SOLID
    from OCP.TopExp import TopExp_Explorer

    shapes = [load(p) for p in paths]
    fused = shapes[0]
    for other in shapes[1:]:
        op = BRepAlgoAPI_Fuse(fused, other)
        op.Build()
        if not op.IsDone():
            print(json.dumps({"done": False}))
            return 1
        fused = op.Shape()
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
    sys.exit(main(sys.argv[1:]))
