// A block from (-10, -10, 0) to (10, 10, 10), united with the paraboloid x² + y² <= 8z,
// z <= 12.5, whose vertex touches the block's bottom face from inside.
let m = model();
let block = m.profile(xy()).rect(20.0, 20.0).sketch("block").extrude(blind(10.0));
let p = m.tool(|m| m.profile(xz()).path(0.0, 0.0).conic_to([5.0, 0.0], [10.0, 12.5], 1.0)
    .line_to(0.0, 12.5).close().sketch("section")
    .revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle()));
block = block.union(p);
m
