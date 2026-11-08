// Two spherical sectors of a ball of radius 10 (cones of half-angle 45 degrees from the centre),
// about z and about z turned 60 degrees about y: their common part.
let m = model();
let first = m.sphere(10.0);
let cone1 = m.tool(|m| m.cone(0.0, 10.0, 10.0));
first = first.intersect(cone1);
let second = m.tool(|m| {
    let s = m.sphere(10.0);
    let c = m.tool(|m| m.cone(0.0, 10.0, 10.0));
    s.intersect(c).rotate(axis_at([0.0, 0.0, 0.0], [0.0, 1.0, 0.0]), 60.0)
});
first = first.intersect(second);
m
