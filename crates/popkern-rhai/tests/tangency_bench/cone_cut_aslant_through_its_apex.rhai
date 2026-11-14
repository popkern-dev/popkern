// The cone of radius 5 at z = 0 with apex (0, 0, 10), less the side of 4x + z = 10 away from the
// axis (the plane through the apex and the line x = 2.5 of the base plane).
let m = model();
let c = m.cone(5.0, 0.0, 10.0);
let n = 1.0 / sqrt(17.0);
let cutter = m.tool(|m| m.profile(plane_at([0.0, 0.0, 10.0], [4.0 * n, 0.0, n])).rect(60.0, 60.0)
    .sketch("cut").extrude(blind(30.0)));
c = c.subtract(cutter);
m
