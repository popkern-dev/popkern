// Two cylinders of radius 5 on the z axis, from z = 0 to 10 and from 6 to 16, united; the second
// turned a quarter turn about z so that the seams of the coincident walls do not meet.
let m = model();
let first = m.profile(xy()).circle(5.0).sketch("first").extrude(blind(10.0));
let second = m.tool(|m| m.profile(xy_at(6.0)).circle(5.0).sketch("second").extrude(blind(10.0))
    .rotate(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), 90.0));
first = first.union(second);
m
