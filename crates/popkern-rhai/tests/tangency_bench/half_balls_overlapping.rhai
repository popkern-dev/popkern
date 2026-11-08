// The half ball z >= 0 (poles on z) united with the half ball x >= 0 (poles on x).
let m = model();
let top = m.sphere(5.0);
let lower = m.tool(|m| m.profile(xy_at(-6.0)).rect(12.0, 12.0).sketch("lower").extrude(blind(6.0)));
top = top.subtract(lower);
let side = m.tool(|m| {
    let b = m.sphere(5.0).rotate(axis_at([0.0, 0.0, 0.0], [0.0, 1.0, 0.0]), 90.0);
    let left = m.tool(|m| m.profile(yz_at(-6.0)).rect(12.0, 12.0).sketch("left").extrude(blind(6.0)));
    b.subtract(left)
});
top = top.union(side);
m
