"""Sew a STEP into a solid like FreeCAD and report its volume."""
import sys
from OCP.STEPControl import STEPControl_Reader
from OCP.IFSelect import IFSelect_RetDone
from OCP.BRepBuilderAPI import BRepBuilderAPI_Sewing, BRepBuilderAPI_MakeSolid
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_SHELL
from OCP.TopoDS import TopoDS
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.BRepCheck import BRepCheck_Analyzer

reader = STEPControl_Reader()
assert reader.ReadFile(sys.argv[1]) == IFSelect_RetDone
reader.TransferRoots()
shape = reader.OneShape()
sew = BRepBuilderAPI_Sewing(1e-6)
sew.Add(shape)
sew.Perform()
sewn = sew.SewedShape()
exp = TopExp_Explorer(sewn, TopAbs_SHELL)
n_shells = 0
volume = None
while exp.More():
    n_shells += 1
    shell = TopoDS.Shell_s(exp.Current())
    mk = BRepBuilderAPI_MakeSolid(shell)
    if mk.IsDone():
        solid = mk.Solid()
        props = GProp_GProps()
        BRepGProp.VolumeProperties_s(solid, props)
        volume = props.Mass()
        ok = BRepCheck_Analyzer(solid).IsValid()
        print(f"shell {n_shells}: solid volume={volume:.1f} valid={ok}")
    exp.Next()
print(f"shells: {n_shells}")
