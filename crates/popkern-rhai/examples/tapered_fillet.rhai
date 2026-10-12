// A variable-radius round: the selector takes both long top edges,
// and each fillet grows from 0.04 to 0.18 along x (a linear law).
let m = model();
let bar = m.profile(xy()).rect(1.5, 0.5).sketch("bar_sketch").extrude(blind(0.4));
bar = bar.fillet(edges().straight().direction("x").at_height(0.4), 0.04,
    #{ law: linear(0.18), label: "taper" });
m.assert_watertight();
m
