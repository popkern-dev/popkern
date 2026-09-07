// Angled fitting: a wedge bracket tapped NORMAL to its inclined face
// (fittings at angles need nothing beyond axis_at — the hole axis
// points into the material, and a screw built on the same axis mates
// by construction), plus a countersunk clearance hole in the base
// and the matching flat-head screw shown alongside (PK-460).
let m = model();

// Wedge: rectangular base, inclined top from z=14 down to z=6.
let bracket = m.profile(xz())
    .path(0.0, 0.0).line_to(24.0, 0.0).line_to(24.0, 6.0).line_to(0.0, 14.0).close()
    .sketch("wedge")
    .extrude(blind(16.0));

// Tap perpendicular to the inclined face (slope (24, -8) in xz →
// inward normal (-8, 0, -24)/25.3), printed-fit clearance.
bracket = bracket.tapped_hole(
    axis_at([12.0, 8.0, 10.0], [-0.316, 0.0, -0.949]),
    "M6",
    blind(6.0),
    #{ clearance: 0.2, label: "angled_tap" });

// A countersunk clearance hole down through the thick end of the base.
bracket = bracket.clearance_hole(
    axis_at([4.0, 8.0, 12.6], [0.0, 0.0, -1.0]),
    "M5",
    through_all(),
    #{ countersink: 0.2, label: "flush_hole" });

// The screw that fits it, standing next to the bracket.
let screw = m.bolt("M5", 12.0, #{ head: "countersunk", thread_length: 7.0, label: "screw" });
screw = screw.translate(32.0, 8.0, 0.0);

m.assert_watertight();
m
