fn ball(m, r) {
    m.sphere(r)
}

// A ball of radius 15 on the origin, united with a cone on the z axis, apex (0, 0, -25), base the
// disc of radius 12 at z = -9. The cone's side is tangent to the
// sphere along the base's rim.
let m = model();
let ice_cream = ball(m, 15.0);
let cone = m.tool(|m| m.cone(0.0, 12.0, 16.0, #{ at: [0.0, 0.0, -25.0] }));
ice_cream = ice_cream.union(cone);
m
