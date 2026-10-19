#!/usr/bin/env python3
"""OpenCASCADE fillet oracle: OCC's own constant-radius fillet of the
curved (non-straight) edges lying at one height, and the result's
validity and volume (PK-888).

Usage:
    <ocp-venv>/bin/python tools/occ-oracle/occ_fillet.py FILE.step Z RADIUS
    <ocp-venv>/bin/python tools/occ-oracle/occ_fillet.py [FILE.step] --angled-hole AX AY AZ Z RADIUS

FILE.step holds the UNFILLETED body. With --angled-hole the body is
built natively in OCC instead (a FILE.step before the flag is ignored), independent of any popkern export: a
40 x 40 x 20 block centred on the z axis with an r 5 cylinder through
(0, 0, 30) along (AX, AY, AZ) subtracted (the PK-888 fixture). Every edge whose curve is not a
line and whose samples all lie within 1e-6 of z = Z is filleted with
RADIUS (an angled hole's top rim: a plane x oblique cylinder ellipse).
The volume is BRepGProp's adaptive integral (1e-9 relative).

Output: one JSON object on stdout:
    edges   how many edges were filleted
    valid   BRepCheck_Analyzer's verdict on the filleted solid
    volume  its volume (mm^3), null when the fillet failed
    ball_worst, fillet_area  (--angled-hole only) how far OCC's fillet
            faces stray from the ideal rolling ball (the worst of
            |d - r| from the ball centre one radius along the normal to
            the plane z = Z and to the hole wall) and their area. OCC
            approximates a non-circular fillet (measured 9.8e-4 at an
            8 deg tilt, 8.2e-3 at 25 deg), so its volume is a reference
            only to within about ball_worst x fillet_area.
Exit code 0 when the file loaded.
"""

import json
import sys


def angled_hole_block(ax, ay, az):
    """The PK-888 block with its angled hole, built natively."""
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder
    from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt

    block = BRepPrimAPI_MakeBox(gp_Pnt(-20.0, -20.0, 0.0), 40.0, 40.0, 20.0).Shape()
    axis = gp_Ax2(gp_Pnt(0.0, 0.0, 30.0), gp_Dir(ax, ay, az))
    tool = BRepPrimAPI_MakeCylinder(axis, 5.0, 60.0).Shape()
    return BRepAlgoAPI_Cut(block, tool).Shape()


def ball_deviation(shape, axis, z, radius):
    """Worst rolling-ball deviation and total area of the fillet faces
    (every face that is not a plane or a cylinder) for the angled-hole
    block: a ball centre one radius along the normal must sit `radius`
    from the plane z = Z and from the r 5 hole wall."""
    import math

    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.BRepGProp import BRepGProp
    from OCP.BRepLProp import BRepLProp_SLProps
    from OCP.BRepTools import BRepTools
    from OCP.GeomAbs import GeomAbs_Cylinder, GeomAbs_Plane
    from OCP.GProp import GProp_GProps
    from OCP.TopAbs import TopAbs_FACE
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS

    norm = math.sqrt(sum(c * c for c in axis))
    a = [c / norm for c in axis]
    worst, area = 0.0, 0.0
    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    while explorer.More():
        face = TopoDS.Face_s(explorer.Current())
        explorer.Next()
        surface = BRepAdaptor_Surface(face)
        if surface.GetType() in (GeomAbs_Plane, GeomAbs_Cylinder):
            continue
        props = GProp_GProps()
        BRepGProp.SurfaceProperties_s(face, props, 1e-9)
        area += props.Mass()
        u1, u2, v1, v2 = BRepTools.UVBounds_s(face)
        for i in range(200):
            for j in range(1, 8):
                u, v = u1 + (u2 - u1) * i / 200, v1 + (v2 - v1) * j / 8
                local = BRepLProp_SLProps(surface, u, v, 1, 1e-9)
                p, n = local.Value(), local.Normal()
                for sign in (1.0, -1.0):
                    c = [p.X() + sign * radius * n.X(), p.Y() + sign * radius * n.Y(),
                         p.Z() + sign * radius * n.Z()]
                    to_plane = z - c[2]
                    w = [c[0], c[1], c[2] - 30.0]
                    along = sum(wi * ai for wi, ai in zip(w, a))
                    to_wall = math.sqrt(sum((wi - ai * along) ** 2 for wi, ai in zip(w, a))) - 5.0
                    if to_plane > 0.0 and to_wall > 0.0:
                        worst = max(worst, abs(to_plane - radius), abs(to_wall - radius))
    return worst, area


def main(path, z, radius, native=None):
    from OCP.BRepAdaptor import BRepAdaptor_Curve
    from OCP.BRepCheck import BRepCheck_Analyzer
    from OCP.BRepFilletAPI import BRepFilletAPI_MakeFillet
    from OCP.BRepGProp import BRepGProp
    from OCP.GeomAbs import GeomAbs_Line
    from OCP.GProp import GProp_GProps
    from OCP.STEPControl import STEPControl_Reader
    from OCP.TopAbs import TopAbs_EDGE, TopAbs_SOLID
    from OCP.TopExp import TopExp, TopExp_Explorer
    from OCP.TopTools import TopTools_IndexedMapOfShape
    from OCP.TopoDS import TopoDS

    if native is not None:
        shape = angled_hole_block(*native)
    else:
        reader = STEPControl_Reader()
        if reader.ReadFile(path) != 1:
            print(json.dumps({"edges": 0, "valid": False, "volume": None}))
            return 1
        reader.TransferRoots()
        shape = reader.OneShape()
    solid = TopExp_Explorer(shape, TopAbs_SOLID).Current()
    fillet = BRepFilletAPI_MakeFillet(solid)
    count = 0
    # One entry per distinct edge (the explorer visits each edge once
    # per face that uses it; Python-side hashes do not identify shapes).
    edges = TopTools_IndexedMapOfShape()
    TopExp.MapShapes_s(solid, TopAbs_EDGE, edges)
    for index in range(1, edges.Extent() + 1):
        edge = TopoDS.Edge_s(edges.FindKey(index))
        curve = BRepAdaptor_Curve(edge)
        if curve.GetType() == GeomAbs_Line:
            continue
        t0, t1 = curve.FirstParameter(), curve.LastParameter()
        if all(
            abs(curve.Value(t0 + (t1 - t0) * k / 16).Z() - z) <= 1e-6 for k in range(17)
        ):
            fillet.Add(radius, edge)
            count += 1
    result = {"edges": count, "valid": False, "volume": None}
    try:
        fillet.Build()
    except Exception:  # OCC raises StdFail_NotDone on a failed fillet
        print(json.dumps(result))
        return 0
    if fillet.IsDone():
        shape = fillet.Shape()
        props = GProp_GProps()
        BRepGProp.VolumeProperties_s(shape, props, 1e-9)
        result["valid"] = BRepCheck_Analyzer(shape).IsValid()
        result["volume"] = abs(props.Mass())
        if native is not None:
            result["ball_worst"], result["fillet_area"] = ball_deviation(
                shape, native, z, radius
            )
    print(json.dumps(result))
    return 0


if __name__ == "__main__":
    if "--angled-hole" in sys.argv:
        at = sys.argv.index("--angled-hole")
        ax, ay, az, z, radius = (float(a) for a in sys.argv[at + 1 : at + 6])
        sys.exit(main(None, z, radius, native=(ax, ay, az)))
    sys.exit(main(sys.argv[1], float(sys.argv[2]), float(sys.argv[3])))
