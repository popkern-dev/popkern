// A knob: a stepped profile revolved about Z, its circular rims
// rounded. Every face is an exact plane, cylinder or cone, so the
// mesh is watertight without repair.
let m = model();
let knob = m.profile(xz())
    .path(0.0, 0.0).line_to(0.45, 0.0).line_to(0.45, 0.25).line_to(0.30, 0.25)
    .line_to(0.30, 0.55).line_to(0.0, 0.55).close()
    .sketch("knob_section")
    .revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle());
knob.fillet(edges().circular().at_height(0.25).or(edges().circular().at_height(0.55)), 0.04,
    #{ label: "rounds" });
m.assert_watertight();
m
