// A ball of radius 5, its poles turned onto (3, 0, 4) and (-3, 0, -4), bored through along z by a
// hole of radius 3: the rims at z = ±4 run through the poles.
let m = model();
let b = m.sphere(5.0).rotate(axis_at([0.0, 0.0, 0.0], [0.0, 1.0, 0.0]), 36.86989764584402);
let bore = m.tool(|m| m.profile(xy_at(-6.0)).circle(3.0).sketch("bore").extrude(blind(12.0)));
b = b.subtract(bore);
m
