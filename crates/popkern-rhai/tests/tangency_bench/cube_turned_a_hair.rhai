// A cube from (-5, -5, 0) to (5, 5, 10), united with a copy turned 0.001 degrees about z.
let m = model();
let cube = m.profile(xy()).rect(10.0, 10.0).sketch("cube").extrude(blind(10.0));
let turned = m.tool(|m| m.profile(xy()).rect(10.0, 10.0).sketch("turned").extrude(blind(10.0))
    .rotate(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), 0.001));
cube = cube.union(turned);
m
