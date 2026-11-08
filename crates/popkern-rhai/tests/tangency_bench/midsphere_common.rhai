fn ball(m, r) {
    m.sphere(r)
}

// A cube of side 10 on the origin, and the ball on the origin of radius 5√2 that touches its
// twelve edges: their common part.
let m = model();
let cube = m.profile(xy_at(-5.0)).rect(10.0, 10.0).sketch("cube").extrude(blind(10.0));
let midsphere = m.tool(|m| ball(m, 5.0 * sqrt(2.0)));
cube = cube.intersect(midsphere);
m
