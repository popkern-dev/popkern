// The cup x² + y² <= 10z, z <= 5 (the parabola x² = 10z, a conic of weight 1 with shoulder
// (2.5, 0), turned about z), united with a ball of radius 5 at (0, 0, 5) whose bottom pole
// touches the cup's vertex to high order.
let m = model();
let cup = m.profile(xz()).path(0.0, 0.0).conic_to([2.5 * sqrt(2.0), 0.0], [5.0 * sqrt(2.0), 5.0], 1.0)
    .line_to(0.0, 5.0).close().sketch("cup")
    .revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle());
let ball = m.tool(|m| m.sphere(5.0).translate(0.0, 0.0, 5.0));
cup = cup.union(ball);
m
