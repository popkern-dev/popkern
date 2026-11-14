// Halves of the torus about z (tube middle 10, tube 3): z >= 0, and x >= 0 of a copy turned a
// quarter turn about z; their union.
let m = model();
let upper = m.profile(xz()).circle(3.0).at(10.0, 0.0).sketch("tube")
    .revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle());
let below = m.tool(|m| m.profile(xy_at(-5.0)).rect(30.0, 30.0).sketch("below").extrude(blind(5.0)));
upper = upper.subtract(below);
let right = m.tool(|m| {
    let t = m.profile(xz()).circle(3.0).at(10.0, 0.0).sketch("tube2")
        .revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle())
        .rotate(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), 90.0);
    let left = m.tool(|m| m.profile(yz_at(-15.0)).rect(30.0, 30.0).sketch("left").extrude(blind(15.0)));
    t.subtract(left)
});
upper = upper.union(right);
m
