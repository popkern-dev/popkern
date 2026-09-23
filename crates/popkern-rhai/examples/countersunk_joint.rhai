// Countersunk joint: a plate with a flush clearance seat and the
// flat-head screw that fits it, side by side. The seat is a plain
// clearance hole with an exact-cone countersink — no thread to cut —
// so the screw's single M5 thread is the only coil in the build.
let m = model();

let plate = m.profile(xy())
    .rect(20.0, 14.0)
    .sketch("plate")
    .extrude(blind(4.0));

plate = plate.clearance_hole(
    axis_at([0.0, 0.0, 4.0], [0.0, 0.0, -1.0]),
    "M5",
    through_all(),
    #{ countersink: 0.2, label: "seat" });

let screw = m.bolt("M5", 10.0, #{
    head: "countersunk",
    thread_length: 5.0,
    tip_chamfer: false,
    label: "screw" });
screw = screw.translate(14.0, 0.0, 0.0);

m.assert_watertight();
m
