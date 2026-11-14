// A ball of radius 10 on the origin and a cylinder of radius 5 on the axis x = 5, y = 0, from
// z = -12 to 12, which touches the ball's sphere from inside at (10, 0, 0): the two united.
// The ball is a half disc revolved about the y axis, its poles at (0, ±10, 0)
// and its seam through (-10, 0, 0).
let m = model();
let ball = m.profile(xy()).path(0.0, -10.0).arc_to([-10.0, 0.0], [0.0, 10.0]).line_to(0.0, -10.0)
    .close().sketch("half_disc").revolve(axis_at([0.0, 0.0, 0.0], [0.0, 1.0, 0.0]), full_circle());
let rod = m.tool(|m| m.profile(xy_at(-12.0)).circle(5.0).at(5.0, 0.0).sketch("rod")
    .extrude(blind(24.0)));
ball = ball.union(rod);
m
