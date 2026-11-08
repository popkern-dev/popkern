// The paraboloid x² + y² <= 8z, z <= 8, and the half-space below z = 1e-4: the tip.
let m = model();
let p = m.profile(xz()).path(0.0, 0.0).conic_to([4.0, 0.0], [8.0, 8.0], 1.0).line_to(0.0, 8.0)
    .close().sketch("section").revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle());
let below = m.tool(|m| m.profile(xy_at(-1.0)).rect(30.0, 30.0).sketch("below").extrude(blind(1.0001)));
p = p.intersect(below);
m
