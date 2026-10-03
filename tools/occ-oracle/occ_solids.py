#!/usr/bin/env python3
"""OpenCASCADE solid oracle: every solid OCC reads from a STEP file,
with its face count and volume (PK-649).

Usage:
    <ocp-venv>/bin/python tools/occ-oracle/occ_solids.py FILE.step

Output: one JSON object on stdout:
    solids  one row per solid (placed occurrences repeat): faces, volume
            (mm^3; the reader converts to millimetres as popkern does)
Exit code 0 when the file loaded.
"""

import json
import sys


def main(path):
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    from OCP.STEPControl import STEPControl_Reader
    from OCP.TopAbs import TopAbs_FACE, TopAbs_SOLID
    from OCP.TopExp import TopExp_Explorer

    reader = STEPControl_Reader()
    if reader.ReadFile(path) != 1:
        print(json.dumps({"solids": []}))
        return 1
    reader.TransferRoots()
    rows = []
    for i in range(1, reader.NbShapes() + 1):
        solids = TopExp_Explorer(reader.Shape(i), TopAbs_SOLID)
        while solids.More():
            solid = solids.Current()
            props = GProp_GProps()
            BRepGProp.VolumeProperties_s(solid, props)
            faces = 0
            explorer = TopExp_Explorer(solid, TopAbs_FACE)
            while explorer.More():
                faces += 1
                explorer.Next()
            rows.append({"faces": faces, "volume": abs(props.Mass())})
            solids.Next()
    print(json.dumps({"solids": rows}))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
