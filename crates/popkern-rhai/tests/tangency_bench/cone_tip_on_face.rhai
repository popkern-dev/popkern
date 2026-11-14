// A block from (-10, -10, 0) to (10, 10, 10), united with a cone on the z axis whose apex is the
// origin and whose base is the disc of radius 5 at z = 15.
let m = model();
let block = m.profile(xy()).rect(20.0, 20.0).sketch("block").extrude(blind(10.0));
let cone = m.tool(|m| m.cone(0.0, 5.0, 15.0));
block = block.union(cone);
m
