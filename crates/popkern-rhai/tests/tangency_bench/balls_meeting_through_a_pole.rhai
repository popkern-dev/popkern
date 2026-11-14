// Balls of radius 5 at the origin and at (4.8, 0, 3.6), poles on the verticals through their
// centres: the circle where they meet passes through a pole of each.
let m = model();
let first = m.sphere(5.0);
let second = m.tool(|m| m.sphere(5.0).translate(4.8, 0.0, 3.6));
first = first.union(second);
m
