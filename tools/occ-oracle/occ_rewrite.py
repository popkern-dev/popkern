#!/usr/bin/env python3
"""OpenCASCADE STEP read-then-write: hands popkern the geometry OCC's
reader produced, so popkern's certifier can audit it (PK-506 round-trip
continuity leg: does continuity survive a foreign reader?).

Usage:
    <ocp-venv>/bin/python tools/occ-oracle/occ_rewrite.py IN.step OUT.step [IN2 OUT2 ...]

Each IN is read with STEPControl_Reader (all roots transferred) and
written with STEPControl_Writer as an AS_IS shape. Output: one JSON
object per pair on stdout with the reader and writer status. Exit
code 0 iff every pair was read and written.
"""

import json
import sys


def rewrite(src, dst):
    from OCP.IFSelect import IFSelect_RetDone
    from OCP.STEPControl import STEPControl_AsIs, STEPControl_Reader, STEPControl_Writer

    reader = STEPControl_Reader()
    status = reader.ReadFile(src)
    if status != IFSelect_RetDone:
        return {"in": src, "ok": False, "stage": "read", "status": str(status)}
    reader.TransferRoots()
    shape = reader.OneShape()
    writer = STEPControl_Writer()
    transfer = writer.Transfer(shape, STEPControl_AsIs)
    if transfer != IFSelect_RetDone:
        return {"in": src, "ok": False, "stage": "transfer", "status": str(transfer)}
    written = writer.Write(dst)
    ok = written == IFSelect_RetDone
    return {"in": src, "out": dst, "ok": ok, "stage": "write", "status": str(written)}


def main(argv):
    pairs = argv[1:]
    if not pairs or len(pairs) % 2:
        print(__doc__, file=sys.stderr)
        return 2
    all_ok = True
    for src, dst in zip(pairs[0::2], pairs[1::2]):
        try:
            result = rewrite(src, dst)
        except Exception as error:  # the caller records the failure
            result = {"in": src, "ok": False, "stage": "exception", "status": repr(error)}
        all_ok &= result["ok"]
        print(json.dumps(result))
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
