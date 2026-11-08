// A ball of radius 5 (poles on z) and the half-space above z = 4.9999: the cap 1e-4 high.
let m = model();
let b = m.sphere(5.0);
let above = m.tool(|m| m.profile(xy_at(4.9999)).rect(20.0, 20.0).sketch("above").extrude(blind(5.0)));
b = b.intersect(above);
m
