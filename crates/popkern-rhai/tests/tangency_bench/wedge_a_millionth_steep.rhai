// A cube from (-5, -5, -10) to (5, 5, 0) and the half-space above the plane z = y·1e-6.
let m = model();
let cube = m.profile(xy_at(-10.0)).rect(10.0, 10.0).sketch("cube").extrude(blind(10.0));
let n = 1.0 / sqrt(1.0 + 1e-12);
let above = m.tool(|m| m.profile(plane_at([0.0, 0.0, 0.0], [0.0, -1e-6 * n, n])).rect(40.0, 40.0)
    .sketch("above").extrude(blind(20.0)));
cube = cube.intersect(above);
m
