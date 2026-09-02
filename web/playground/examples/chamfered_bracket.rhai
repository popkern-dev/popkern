// An L-bracket with three chamfer kinds: equal on the outer corner,
// two-distance on the inner corner, distance-and-angle on the ends.
let m = model();
let bracket = m.profile(yz())
    .path(0.0, 0.0).line_to(1.0, 0.0).line_to(1.0, 0.15).line_to(0.15, 0.15)
    .line_to(0.15, 0.8).line_to(0.0, 0.8).close()
    .sketch("section")
    .extrude(blind(0.6));
bracket = bracket.chamfer(edges().straight().direction("x").at_height(0.0), 0.05,
    #{ label: "outer" });
bracket = bracket.chamfer(edges().straight().direction("x").at_height(0.8), 0.04,
    #{ kind: distance_angle(faces().facing("z").at_height(0.8), 30.0), label: "top" });
bracket = bracket.chamfer(edges().straight().direction("x").at_height(0.15).allow_empty(), 0.03,
    #{ kind: two_distance(faces().facing("z").at_height(0.15), 0.06), label: "step" });
m.assert_watertight();
m
