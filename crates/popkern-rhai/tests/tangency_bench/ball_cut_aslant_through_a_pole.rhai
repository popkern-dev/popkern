// A ball of radius 5 (poles on z) less the cap beyond the plane x + z = 5, whose rim passes
// through the pole (0, 0, 5).
let m = model();
let b = m.sphere(5.0);
let s = 1.0 / sqrt(2.0);
let cutter = m.tool(|m| m.profile(plane_at([2.5, 0.0, 2.5], [s, 0.0, s])).rect(40.0, 40.0)
    .sketch("cut").extrude(blind(20.0)));
b = b.subtract(cutter);
m
