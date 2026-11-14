// A torus about z (centre circle 10, tube 3), united with a drum of radius 13 on the z axis from
// z = -1 to 1, which it touches from inside along its outer equator.
let m = model();
let torus = m.profile(xz()).circle(3.0).at(10.0, 0.0).sketch("tube")
    .revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle());
let drum = m.tool(|m| m.profile(xy_at(-1.0)).circle(13.0).sketch("drum").extrude(blind(2.0)));
torus = torus.union(drum);
m
