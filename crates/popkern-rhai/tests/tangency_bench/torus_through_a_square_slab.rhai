// A torus about z (tube middle 10, tube 3, seam turned 45 degrees) united with the slab
// (-13, -13, -1) to (13, 13, 1), whose sides it touches at their middles.
let m = model();
let torus = m.profile(xz()).circle(3.0).at(10.0, 0.0).sketch("tube")
    .revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle())
    .rotate(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), 45.0);
let slab = m.tool(|m| m.profile(xy_at(-1.0)).rect(26.0, 26.0).sketch("slab").extrude(blind(2.0)));
torus = torus.union(slab);
m
