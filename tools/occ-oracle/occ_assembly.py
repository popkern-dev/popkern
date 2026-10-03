#!/usr/bin/env python3
"""OpenCASCADE assembly oracle: how an independent reader (XCAF) sees
a STEP product structure (PK-629).

Usage:
    <ocp-venv>/bin/python tools/occ-oracle/occ_assembly.py FILE.step

Output: one JSON object on stdout:
    reader_status   load status string
    check_messages  Part 21 reader check messages (fail/warn texts)
    parts           number of distinct part (leaf) shapes in the shape tool
    instances       one row per placed leaf occurrence, depth-first:
                    path   instance names from the root ("stage/block_1")
                    part   the referred part's product name
                    bbox   located bounding box [[xmin,ymin,zmin],[xmax,ymax,zmax]]
                    solids solid count of the located shape
Exit code 0 when the file loaded.
"""

import json
import sys


def label_name(label):
    from OCP.TDataStd import TDataStd_Name
    from OCP.TCollection import TCollection_ExtendedString

    attr = TDataStd_Name()
    if label.FindAttribute(TDataStd_Name.GetID_s(), attr):
        text = attr.Get()
        return "".join(text.Value(i) for i in range(1, text.Length() + 1))
    return ""


def bbox(shape):
    from OCP.Bnd import Bnd_Box
    from OCP.BRepBndLib import BRepBndLib

    box = Bnd_Box()
    box.SetGap(0.0)
    BRepBndLib.AddOptimal_s(shape, box, False, False)
    lo = box.CornerMin()
    hi = box.CornerMax()
    return [[lo.X(), lo.Y(), lo.Z()], [hi.X(), hi.Y(), hi.Z()]]


def solids(shape):
    from OCP.TopAbs import TopAbs_SOLID
    from OCP.TopExp import TopExp_Explorer

    n = 0
    exp = TopExp_Explorer(shape, TopAbs_SOLID)
    while exp.More():
        n += 1
        exp.Next()
    return n


def main(path):
    from OCP.IFSelect import IFSelect_RetDone
    from OCP.STEPCAFControl import STEPCAFControl_Reader
    from OCP.TCollection import TCollection_ExtendedString
    from OCP.TDF import TDF_Label, TDF_LabelSequence
    from OCP.TDocStd import TDocStd_Document
    from OCP.TopLoc import TopLoc_Location
    from OCP.XCAFDoc import XCAFDoc_DocumentTool

    doc = TDocStd_Document(TCollection_ExtendedString("pk-assembly"))
    reader = STEPCAFControl_Reader()
    reader.SetNameMode(True)
    status = reader.ReadFile(path)
    out = {"reader_status": str(status), "check_messages": [], "parts": 0, "instances": []}
    if status != IFSelect_RetDone:
        print(json.dumps(out))
        return 1
    reader.Transfer(doc)
    tool = XCAFDoc_DocumentTool.ShapeTool_s(doc.Main())

    parts = set()

    def walk(label, location, path):
        if tool.IsAssembly_s(label):
            comps = TDF_LabelSequence()
            tool.GetComponents_s(label, comps, False)
            for i in range(1, comps.Length() + 1):
                comp = comps.Value(i)
                referred = TDF_Label()
                tool.GetReferredShape_s(comp, referred)
                loc = location.Multiplied(tool.GetLocation_s(comp))
                walk(referred, loc, path + [label_name(comp) or label_name(referred)])
        else:
            parts.add(label.Tag())
            shape = tool.GetShape_s(label).Located(location)
            out["instances"].append(
                {
                    "path": "/".join(path),
                    "part": label_name(label),
                    "bbox": bbox(shape),
                    "solids": solids(shape),
                }
            )

    roots = TDF_LabelSequence()
    tool.GetFreeShapes(roots)
    for i in range(1, roots.Length() + 1):
        walk(roots.Value(i), TopLoc_Location(), [])
    out["parts"] = len(parts)
    print(json.dumps(out))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
