// A cone, base the disc of radius 5 at z = 0, apex (0, 0, 10); everything with x < 0 cut away.
let m = model();
let c = m.cone(5.0, 0.0, 10.0);
let left = m.tool(|m| m.profile(yz_at(-10.0)).rect(30.0, 30.0).sketch("left").extrude(blind(10.0)));
c = c.subtract(left);
m
