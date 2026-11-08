fn ball(m, r) {
    m.sphere(r)
}

// A cube of side 10 from the origin, united with its inscribed ball, radius 5 at (5, 5, 5).
let m = model();
let cube = m.profile(xy()).rect(10.0, 10.0).at(5.0, 5.0).sketch("cube").extrude(blind(10.0));
let inscribed = m.tool(|m| ball(m, 5.0).translate(5.0, 5.0, 5.0));
cube = cube.union(inscribed);
m
