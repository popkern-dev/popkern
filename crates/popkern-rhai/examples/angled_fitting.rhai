// Angled fitting: a wedge bracket tapped NORMAL to its inclined face
// (fittings at angles need nothing beyond axis_at — the hole axis
// points into the material, and a screw built on the same axis mates
// by construction), plus a countersunk clearance hole in the base
// and the matching flat-head screw shown alongside (PK-460).
let m = model();

// Wedge: rectangular base, inclined top from z=14 down to z=6. The xz
// sketch plane's normal is -y, so the wedge spans y from -16 to 0 and
// its mid-width is y = -8 (the holes sat at y = +8, outside the part,
// and cut nothing while serving).
let bracket = m.profile(xz())
    .path(0.0, 0.0).line_to(24.0, 0.0).line_to(24.0, 6.0).line_to(0.0, 14.0).close()
    .sketch("wedge")
    .extrude(blind(16.0));

// Tap perpendicular to the inclined face (slope (24, -8) in xz →
// inward normal (-8, 0, -24)/|…|, computed exactly: the rounded
// (-0.316, 0, -0.949) sits 0.013° off the normal), printed-fit clearance.
let n_len = sqrt(8.0 * 8.0 + 24.0 * 24.0);
bracket = bracket.tapped_hole(
    axis_at([12.0, -8.0, 10.0], [-8.0 / n_len, 0.0, -24.0 / n_len]),
    "M6",
    blind(6.0),
    #{ clearance: 0.2, label: "angled_tap" });

// A countersunk clearance hole up through the thick end of the base,
// countersunk in the flat underside (a countersink wants a flat face:
// on the sloped top at x = 4 the ~10.4 mm cone bit the end edge).
bracket = bracket.clearance_hole(
    axis_at([6.0, -8.0, 0.0], [0.0, 0.0, 1.0]),
    "M5",
    through_all(),
    #{ countersink: 0.2, label: "flush_hole" });

// The screw that fits it, standing next to the bracket.
let screw = m.bolt("M5", 12.0, #{ head: "countersunk", thread_length: 7.0, label: "screw" });
screw = screw.translate(32.0, -8.0, 0.0);

m.assert_watertight();
m
