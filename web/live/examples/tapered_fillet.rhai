// A variable-radius round: the fillet grows from 0.04 to 0.18 along
// the top front edge (a linear law), a constant round takes the back.
let m = model();
let bar = m.profile(xy()).rect(1.5, 0.5).sketch("bar_sketch").extrude(blind(0.4));
bar = bar.fillet(edges().straight().direction("x").at_height(0.4), 0.04,
    #{ law: linear(0.18), label: "taper" });
m.assert_watertight();
m
