"""OpenCASCADE oracle: import a STEP like FreeCAD and report per-face trim health."""
import sys
from OCP.STEPControl import STEPControl_Reader
from OCP.IFSelect import IFSelect_RetDone
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_FACE, TopAbs_EDGE, TopAbs_WIRE
from OCP.TopoDS import TopoDS
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.GeomAbs import GeomAbs_Plane, GeomAbs_Cylinder, GeomAbs_Sphere, GeomAbs_Torus
from OCP.ShapeFix import ShapeFix_Shape
from OCP.BRepMesh import BRepMesh_IncrementalMesh
from OCP.ShapeAnalysis import ShapeAnalysis_FreeBounds

path = sys.argv[1]
reader = STEPControl_Reader()
status = reader.ReadFile(path)
assert status == IFSelect_RetDone, "read failed"
reader.TransferRoots()
shape = reader.OneShape()

kind_names = {GeomAbs_Plane: "Plane", GeomAbs_Cylinder: "Cylinder",
              GeomAbs_Sphere: "Sphere", GeomAbs_Torus: "Torus"}

print("== overall validity:", BRepCheck_Analyzer(shape).IsValid())

exp = TopExp_Explorer(shape, TopAbs_FACE)
i = 0
while exp.More():
    face = TopoDS.Face_s(exp.Current())
    surf = BRepAdaptor_Surface(face)
    kind = kind_names.get(surf.GetType(), str(surf.GetType()))
    props = GProp_GProps()
    BRepGProp.SurfaceProperties_s(face, props)
    area = props.Mass()
    valid = BRepCheck_Analyzer(face).IsValid()
    # count edges
    ne = 0
    ee = TopExp_Explorer(face, TopAbs_EDGE)
    while ee.More():
        ne += 1
        ee.Next()
    extra = ""
    if kind == "Cylinder":
        c = surf.Cylinder()
        loc = c.Location()
        extra = f" r={c.Radius():.4f} at ({loc.X():.3f},{loc.Y():.3f},{loc.Z():.3f})"
    if kind == "Torus":
        t = surf.Torus()
        extra = f" R={t.MajorRadius():.4f} r={t.MinorRadius():.4f}"
    print(f"face[{i}] {kind:8s} area={area:9.5f} edges={ne} valid={valid}{extra}")
    i += 1
    exp.Next()

# free boundary census after meshing (the visual-gap oracle)
BRepMesh_IncrementalMesh(shape, 0.005, False, 0.2, False)
fb = ShapeAnalysis_FreeBounds(shape)
closed = fb.GetClosedWires()
open_w = fb.GetOpenWires()
def count_wires(comp):
    n = 0
    e = TopExp_Explorer(comp, TopAbs_WIRE)
    while e.More():
        n += 1
        e.Next()
    return n
print("== free bounds: closed wires:", count_wires(closed), " open wires:", count_wires(open_w))

# Locate free-bound wires (appended): print each wire's bbox.
from OCP.Bnd import Bnd_Box
from OCP.BRepBndLib import BRepBndLib
for label, compound in (("closed", fb.GetClosedWires()), ("open", fb.GetOpenWires())):
    exp = TopExp_Explorer(compound, TopAbs_EDGE)
    while exp.More():
        box = Bnd_Box()
        BRepBndLib.Add_s(exp.Current(), box)
        if not box.IsVoid():
            mn = box.CornerMin()
            mx = box.CornerMax()
            print(f"  {label} free edge bbox z=[{mn.Z():.2f},{mx.Z():.2f}] x=[{mn.X():.2f},{mx.X():.2f}]")
        exp.Next()
