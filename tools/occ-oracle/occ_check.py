#!/usr/bin/env python3
"""OpenCASCADE STEP oracle: reader checks, BRepCheck census, free bounds.

The FreeCAD-equivalent acceptance gate for popkern STEP exports
(see docs: popkern round-trip audits pass defects that live in
STEP-consumer conventions; this oracle surfaces them).

Usage:
    <ocp-venv>/bin/python tools/occ-oracle/occ_check.py [--self-interference] FILE.step [FILE2.step ...]
    # venv recreate: python3 -m venv ~/venvs/ocp && ~/venvs/ocp/bin/pip install cadquery-ocp

Output: one JSON object per file on stdout with:
    reader_status        load status string
    reader_check_messages  Part 21 syntax/attribute defects (f3d silently heals these)
    n_faces / n_shells / surface_kinds   census
    overall_valid        BRepCheck_Analyzer on the whole shape
    face_failures        per-face BRepCheck failure names + surface kind + area
    free_bound_wires     closed/open free-bound wire counts after meshing
    self_interference    (with --self-interference) BRepAlgoAPI_Check
                         self-interference + small-edge verdict and fault
                         statuses — the boolean pre-check a consumer runs
                         before any operation; it catches import-side
                         tolerance inflation that BRepCheck accepts (PK-524)
Exit code 0 iff every file passes: valid, no reader checks, 0 open free bounds.
"""

import json
import sys


def check_file(path, self_interference=False):
    from OCP.BRep import BRep_Tool
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.BRepCheck import BRepCheck_Analyzer, BRepCheck_ListOfStatus
    from OCP.BRepGProp import BRepGProp
    from OCP.BRepMesh import BRepMesh_IncrementalMesh
    from OCP.GProp import GProp_GProps
    from OCP.GeomAbs import GeomAbs_SurfaceType
    from OCP.IFSelect import IFSelect_ItemsByEntity, IFSelect_RetDone
    from OCP.Interface import Interface_CheckIterator
    from OCP.STEPControl import STEPControl_Reader
    from OCP.ShapeAnalysis import ShapeAnalysis_FreeBounds
    from OCP.TopAbs import TopAbs_FACE, TopAbs_SHELL, TopAbs_SOLID, TopAbs_WIRE
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS

    report = {"file": path}

    reader = STEPControl_Reader()
    status = reader.ReadFile(path)
    report["reader_status"] = str(status)
    if status != IFSelect_RetDone:
        report["fatal"] = "read failed"
        return report, False

    # Part 21 conformance messages the reader logged while loading.
    check_messages = []
    checks: Interface_CheckIterator = reader.WS().ModelCheckList()
    checks.Start()
    while checks.More():
        check = checks.Value()
        for i in range(1, check.NbFails() + 1):
            check_messages.append("FAIL: " + check.CFail(i))
        for i in range(1, check.NbWarnings() + 1):
            check_messages.append("WARN: " + check.CWarning(i))
        checks.Next()
    report["reader_check_messages"] = check_messages

    reader.TransferRoots()
    shape = reader.OneShape()

    surface_kind_names = {
        int(GeomAbs_SurfaceType.GeomAbs_Plane): "Plane",
        int(GeomAbs_SurfaceType.GeomAbs_Cylinder): "Cylinder",
        int(GeomAbs_SurfaceType.GeomAbs_Cone): "Cone",
        int(GeomAbs_SurfaceType.GeomAbs_Sphere): "Sphere",
        int(GeomAbs_SurfaceType.GeomAbs_Torus): "Torus",
        int(GeomAbs_SurfaceType.GeomAbs_BezierSurface): "Bezier",
        int(GeomAbs_SurfaceType.GeomAbs_BSplineSurface): "BSpline",
        int(GeomAbs_SurfaceType.GeomAbs_SurfaceOfRevolution): "Revolution",
        int(GeomAbs_SurfaceType.GeomAbs_SurfaceOfExtrusion): "Extrusion",
        int(GeomAbs_SurfaceType.GeomAbs_OffsetSurface): "Offset",
        int(GeomAbs_SurfaceType.GeomAbs_OtherSurface): "Other",
    }

    n_faces = 0
    surface_kinds = {}
    face_failures = []
    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    while explorer.More():
        face = TopoDS.Face_s(explorer.Current())
        n_faces += 1
        kind = surface_kind_names.get(
            int(BRepAdaptor_Surface(face).GetType()), "Unknown"
        )
        surface_kinds[kind] = surface_kinds.get(kind, 0) + 1

        analyzer = BRepCheck_Analyzer(face)
        if not analyzer.IsValid():
            props = GProp_GProps()
            BRepGProp.SurfaceProperties_s(face, props)
            failures = set()
            result = analyzer.Result(face)
            if result is not None:
                statuses: BRepCheck_ListOfStatus = result.Status()
                for st in statuses:
                    failures.add(str(st))
            # Sub-shape statuses (wires carry IntersectingWires etc.)
            sub = TopExp_Explorer(face, TopAbs_WIRE)
            while sub.More():
                sub_result = analyzer.Result(sub.Current())
                if sub_result is not None:
                    for st in sub_result.Status():
                        failures.add(str(st))
                sub.Next()
            failures.discard("BRepCheck_Status.BRepCheck_NoError")
            face_failures.append(
                {
                    "face_index": n_faces - 1,
                    "surface": kind,
                    "area": props.Mass(),
                    "failures": sorted(failures),
                }
            )
        explorer.Next()
    report["n_faces"] = n_faces
    report["surface_kinds"] = surface_kinds
    report["face_failures"] = face_failures

    n_shells = 0
    shell_exp = TopExp_Explorer(shape, TopAbs_SHELL)
    while shell_exp.More():
        n_shells += 1
        shell_exp.Next()
    report["n_shells"] = n_shells

    report["overall_valid"] = BRepCheck_Analyzer(shape).IsValid()

    BRepMesh_IncrementalMesh(shape, 0.1, False, 0.5, True)
    free_bounds = ShapeAnalysis_FreeBounds(shape)
    closed = free_bounds.GetClosedWires()
    opened = free_bounds.GetOpenWires()

    def count_wires(compound):
        n = 0
        exp = TopExp_Explorer(compound, TopAbs_WIRE)
        while exp.More():
            n += 1
            exp.Next()
        return n

    report["free_bound_wires"] = {
        "closed": count_wires(closed),
        "open": count_wires(opened),
    }

    si_ok = True
    if self_interference:
        from OCP.BRepAlgoAPI import BRepAlgoAPI_Check

        # Per SOLID: the check's self-interference semantics are per
        # shape, and a multi-body compound would report contact
        # between separate bodies (a layout fact, not a defect).
        solids = []
        solid_exp = TopExp_Explorer(shape, TopAbs_SOLID)
        while solid_exp.More():
            solids.append(solid_exp.Current())
            solid_exp.Next()
        faults = set()
        for body in solids or [shape]:
            check = BRepAlgoAPI_Check(body, True, True)
            si_ok = si_ok and check.IsValid()
            faults |= {str(res.GetCheckStatus()) for res in check.Result()}
        report["self_interference"] = {
            "valid": si_ok,
            "bodies": len(solids) or 1,
            "faults": sorted(faults),
        }

    passed = (
        report["overall_valid"]
        and not check_messages
        and report["free_bound_wires"]["open"] == 0
        and report["free_bound_wires"]["closed"] == 0
        and si_ok
    )
    report["passed"] = passed
    return report, passed


def main():
    if len(sys.argv) < 2:
        print(__doc__, file=sys.stderr)
        return 2
    self_interference = "--self-interference" in sys.argv[1:]
    paths = [a for a in sys.argv[1:] if a != "--self-interference"]
    if not paths:
        print(__doc__, file=sys.stderr)
        return 2
    all_passed = True
    for path in paths:
        try:
            report, passed = check_file(path, self_interference)
        except Exception as exc:  # report per-file, keep sweeping
            report, passed = {"file": path, "fatal": repr(exc)}, False
        all_passed = all_passed and passed
        print(json.dumps(report, indent=2))
    return 0 if all_passed else 1


if __name__ == "__main__":
    sys.exit(main())
