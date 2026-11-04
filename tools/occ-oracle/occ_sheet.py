#!/usr/bin/env python3
"""OpenCASCADE sheet x solid oracle (PK-978): the same box and sheet
popkern builds, built here from their parameters, so the two kernels
share neither boolean code nor an exchange format.

Usage:
    <ocp-venv>/bin/python tools/occ-oracle/occ_sheet.py SPEC_JSON

SPEC_JSON (one argument):
    solid      {"kind": "box", "corners": [x0, y0, z0, x1, y1, z1]}
               {"kind": "sphere", "centre": [x, y, z], "radius": r}
    sheet      {"kind": "plane", "x": [x0, x1], "y": [y0, y1], "z": h}
               {"kind": "tube", "centre": [x, y], "radius": r, "z": [z0, z1]}
               either with an optional "tilt_x": radians, a rotation
               of the sheet alone about the x axis before placement
    placement  {"axis": [ax, ay, az], "angle": radians,
                "translation": [tx, ty, tz], "scale": s}
               applied as scale about the origin, then the rotation,
               then the translation (popkern's order in the test).

Output: one JSON object on stdout:
    done           every operation completed
    common_area    area of Common(solid, sheet): the sheet inside
    cut_area       area of Cut(sheet, solid): the sheet outside
    split_volumes  volumes of Splitter(solid | sheet)'s solids, ascending
"""

import json
import math
import sys


def placed(shape, placement):
    from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
    from OCP.gp import gp_Ax1, gp_Dir, gp_Pnt, gp_Trsf, gp_Vec

    scale = gp_Trsf()
    scale.SetScale(gp_Pnt(0, 0, 0), placement.get("scale", 1.0))
    rotation = gp_Trsf()
    axis = placement.get("axis", [0, 0, 1])
    rotation.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(*axis)), placement.get("angle", 0.0))
    translation = gp_Trsf()
    translation.SetTranslation(gp_Vec(*placement.get("translation", [0, 0, 0])))
    total = translation.Multiplied(rotation).Multiplied(scale)
    return BRepBuilderAPI_Transform(shape, total, True).Shape()


def sheet_shape(spec):
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace
    from OCP.gp import gp_Ax3, gp_Cylinder, gp_Dir, gp_Pln, gp_Pnt

    tilt = {"axis": [1, 0, 0], "angle": spec.get("tilt_x", 0.0)}
    if spec["kind"] == "plane":
        (x0, x1), (y0, y1), z = spec["x"], spec["y"], spec["z"]
        plane = gp_Pln(gp_Pnt(0, 0, z), gp_Dir(0, 0, 1))
        return placed(BRepBuilderAPI_MakeFace(plane, x0, x1, y0, y1).Face(), tilt)
    if spec["kind"] == "tube":
        (cx, cy), r, (z0, z1) = spec["centre"], spec["radius"], spec["z"]
        cylinder = gp_Cylinder(gp_Ax3(gp_Pnt(cx, cy, z0), gp_Dir(0, 0, 1)), r)
        face = BRepBuilderAPI_MakeFace(cylinder, 0.0, 2.0 * math.pi, 0.0, z1 - z0).Face()
        return placed(face, tilt)
    raise SystemExit(json.dumps({"done": False, "error": f"unknown sheet {spec['kind']}"}))


def solid_shape(spec):
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeSphere
    from OCP.gp import gp_Pnt

    if spec["kind"] == "box":
        x0, y0, z0, x1, y1, z1 = spec["corners"]
        return BRepPrimAPI_MakeBox(gp_Pnt(x0, y0, z0), gp_Pnt(x1, y1, z1)).Shape()
    if spec["kind"] == "sphere":
        return BRepPrimAPI_MakeSphere(gp_Pnt(*spec["centre"]), spec["radius"]).Shape()
    raise SystemExit(json.dumps({"done": False, "error": f"unknown solid {spec['kind']}"}))


def area(shape):
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps

    props = GProp_GProps()
    BRepGProp.SurfaceProperties_s(shape, props, 1e-9)
    return props.Mass()


def volumes(shape):
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    from OCP.TopAbs import TopAbs_SOLID
    from OCP.TopExp import TopExp_Explorer

    out = []
    explorer = TopExp_Explorer(shape, TopAbs_SOLID)
    while explorer.More():
        props = GProp_GProps()
        BRepGProp.VolumeProperties_s(explorer.Current(), props, 1e-9)
        out.append(props.Mass())
        explorer.Next()
    return sorted(out)


def main(spec):
    from OCP.BOPAlgo import BOPAlgo_Splitter
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut

    placement = spec.get("placement", {})
    box = placed(solid_shape(spec["solid"]), placement)
    sheet = placed(sheet_shape(spec["sheet"]), placement)

    common = BRepAlgoAPI_Common(box, sheet)
    cut = BRepAlgoAPI_Cut(sheet, box)
    splitter = BOPAlgo_Splitter()
    splitter.AddArgument(box)
    splitter.AddTool(sheet)
    splitter.Perform()
    done = common.IsDone() and cut.IsDone() and not splitter.HasErrors()
    print(
        json.dumps(
            {
                "done": bool(done),
                "common_area": area(common.Shape()) if common.IsDone() else None,
                "cut_area": area(cut.Shape()) if cut.IsDone() else None,
                "split_volumes": volumes(splitter.Shape()) if not splitter.HasErrors() else None,
            }
        )
    )


if __name__ == "__main__":
    main(json.loads(sys.argv[1]))
