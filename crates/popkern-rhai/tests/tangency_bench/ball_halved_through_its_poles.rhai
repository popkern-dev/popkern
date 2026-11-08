// A ball of radius 5 (poles on z), everything with x < 0 cut away by a block.
let m = model();
let b = m.sphere(5.0);
let left = m.tool(|m| m.profile(yz_at(-10.0)).rect(20.0, 20.0).sketch("left").extrude(blind(10.0)));
b = b.subtract(left);
m
