# STEP oracle — OpenCASCADE validation for exported STEP files

FreeCAD-equivalent gates for popkern STEP exports. Watertight
tessellation and a popkern round-trip do NOT imply OpenCASCADE can
sew the shell (measured: cap rim edges exported as forked face-local
copies read as an open pipe shell + detached socket prisms in
FreeCAD, fixed in 2c02bb62). Every shipped STEP should pass BOTH:

- `face_health.py <file.step>` — per-face BRepCheck validity, surface
  census, and free-bounds count (must be 0/0).
- `sew_to_solid.py <file.step>` — BRepBuilderAPI_Sewing →
  MakeSolid → GProp volume. Must yield ONE shell, a VALID solid, and
  a volume matching the analytic expectation (catches inverted trims
  and silent no-ops that face checks miss).

Setup: `python3 -m venv occenv && occenv/bin/pip install cadquery-ocp`
