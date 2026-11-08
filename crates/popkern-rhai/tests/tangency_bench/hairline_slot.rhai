// A cube of side 10 from the origin, with a slot 1e-6 wide, x from 5 to 5.000001, cut down to
// z = 5 by a block that runs past the cube's front, back and top.
let m = model();
let cube = m.profile(xy()).rect(10.0, 10.0).at(5.0, 5.0).sketch("cube").extrude(blind(10.0));
let saw = m.tool(|m| m.profile(xy_at(5.0)).rect(0.000001, 12.0).at(5.0000005, 5.0)
    .sketch("saw").extrude(blind(6.0)));
cube = cube.subtract(saw);
m
