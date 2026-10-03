#!/usr/bin/env python3
"""Write a small multi-part STEP assembly with OpenCASCADE's XCAF writer
(PK-631): a third-party writer's product structure for popkern's
importer to read — nested sub-assembly, shared parts placed several
times, instance and product names.

    <ocp-venv>/bin/python tools/occ-oracle/occ_write_assembly.py OUT.step

Structure (lengths in mm):
    bench                         (root)
      top        : plate 200x80x10    at z = 60
      leg_front  : leg  (cylinder r8 h60) at ( 85,  30, 0)
      leg_back   : leg                  at (-85,  30, 0)
      pair       : leg_pair (sub-assembly) at (0, -30, 0)
        left     : leg at ( 85, 0, 0)
        right    : leg at (-85, 0, 0) rotated 90 deg about z
"""

import sys


def main(out):
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder
    from OCP.IFSelect import IFSelect_RetDone
    from OCP.STEPCAFControl import STEPCAFControl_Writer
    from OCP.STEPControl import STEPControl_AsIs
    from OCP.TCollection import TCollection_ExtendedString
    from OCP.TDataStd import TDataStd_Name
    from OCP.TDocStd import TDocStd_Document
    from OCP.TopLoc import TopLoc_Location
    from OCP.XCAFDoc import XCAFDoc_DocumentTool
    from OCP.gp import gp_Ax1, gp_Dir, gp_Pnt, gp_Trsf, gp_Vec

    doc = TDocStd_Document(TCollection_ExtendedString("bench"))
    tool = XCAFDoc_DocumentTool.ShapeTool_s(doc.Main())

    def name(label, text):
        TDataStd_Name.Set_s(label, TCollection_ExtendedString(text))

    def at(x, y, z, turn_deg=0.0):
        t = gp_Trsf()
        if turn_deg:
            t.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), turn_deg * 3.141592653589793 / 180.0)
        move = gp_Trsf()
        move.SetTranslation(gp_Vec(x, y, z))
        return TopLoc_Location(move.Multiplied(t))

    plate = tool.AddShape(BRepPrimAPI_MakeBox(gp_Pnt(-100, -40, 0), 200, 80, 10).Shape(), False)
    name(plate, "plate")
    leg = tool.AddShape(BRepPrimAPI_MakeCylinder(8, 60).Shape(), False)
    name(leg, "leg")
    pair = tool.NewShape()
    name(pair, "leg_pair")
    name(tool.AddComponent(pair, leg, at(85, 0, 0)), "left")
    name(tool.AddComponent(pair, leg, at(-85, 0, 0, 90.0)), "right")
    root = tool.NewShape()
    name(root, "bench")
    name(tool.AddComponent(root, plate, at(0, 0, 60)), "top")
    name(tool.AddComponent(root, leg, at(85, 30, 0)), "leg_front")
    name(tool.AddComponent(root, leg, at(-85, 30, 0)), "leg_back")
    name(tool.AddComponent(root, pair, at(0, -30, 0)), "pair")
    tool.UpdateAssemblies()

    writer = STEPCAFControl_Writer()
    writer.SetNameMode(True)
    if not writer.Transfer(doc, STEPControl_AsIs):
        return 1
    return 0 if writer.Write(out) == IFSelect_RetDone else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
