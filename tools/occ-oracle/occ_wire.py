#!/usr/bin/env python3
"""OpenCASCADE wire x solid oracle (PK-1068): the same solid and wire
popkern builds, built here from their parameters (no exchange format).

Usage:
    <ocp-venv>/bin/python tools/occ-oracle/occ_wire.py SPEC_JSON

SPEC_JSON (one argument):
    solid      as occ_sheet.py: {"kind": "box", "corners": [...]} or
               {"kind": "sphere", "centre": [...], "radius": r}
    wire       {"kind": "polyline", "points": [[x, y, z], ...]}
               {"kind": "circle", "centre": [x, y, z], "radius": r}
               (a full circle in a z plane)
    placement  as occ_sheet.py (scale, then rotation, then translation)

Output: one JSON object on stdout:
    done           every operation completed
    common_length  length of Common(wire, solid): the wire inside
    cut_length     length of Cut(wire, solid): the wire outside
"""

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from occ_sheet import placed, solid_shape  # noqa: E402


def wire_shape(spec):
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire
    from OCP.gp import gp_Ax2, gp_Circ, gp_Dir, gp_Pnt

    if spec["kind"] == "polyline":
        wire = BRepBuilderAPI_MakeWire()
        pts = [gp_Pnt(*p) for p in spec["points"]]
        for a, b in zip(pts, pts[1:]):
            wire.Add(BRepBuilderAPI_MakeEdge(a, b).Edge())
        return wire.Wire()
    if spec["kind"] == "circle":
        circle = gp_Circ(gp_Ax2(gp_Pnt(*spec["centre"]), gp_Dir(0, 0, 1)), spec["radius"])
        return BRepBuilderAPI_MakeEdge(circle).Edge()
    raise SystemExit(json.dumps({"done": False, "error": f"unknown wire {spec['kind']}"}))


def length(shape):
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps

    props = GProp_GProps()
    BRepGProp.LinearProperties_s(shape, props)
    return props.Mass()


def main(spec):
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut

    placement = spec.get("placement", {})
    solid = placed(solid_shape(spec["solid"]), placement)
    wire = placed(wire_shape(spec["wire"]), placement)
    common = BRepAlgoAPI_Common(wire, solid)
    cut = BRepAlgoAPI_Cut(wire, solid)
    done = common.IsDone() and cut.IsDone()
    print(
        json.dumps(
            {
                "done": bool(done),
                "common_length": length(common.Shape()) if common.IsDone() else None,
                "cut_length": length(cut.Shape()) if cut.IsDone() else None,
            }
        )
    )


if __name__ == "__main__":
    main(json.loads(sys.argv[1]))
