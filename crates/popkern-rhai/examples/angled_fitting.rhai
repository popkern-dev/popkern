// Angled fitting, as an ASSEMBLY: a wedge bracket that holds a part
// at an angle. The M6 tap is cut NORMAL to the inclined face, and the
// cap screw that fills it stands on the very same axis — a hole axis
// points into the material and a fastener's insertion frame runs +Z
// down its shank from the bearing plane (PK-460), so placing one on
// the other lines the threads up. The screw bears on a washer: its
// last turn under the head runs out to the full major diameter, which
// a head seated straight on the tap's face would drive into the first
// thread (PK-849). The flat flange at the thin end is screwed down to the
// bench through two countersunk clearance holes, flat-head screws
// flush in their seats.
let a = assembly("angled_fitting");

// Wedge: inclined top from z = 14 at x = 0 down to z = 5 at x = 24,
// then a 5 mm flange out to x = 44. The xz sketch plane's normal is
// -y, so the bracket spans y from -24 to 0; its mid-width is y = -12.
let bracket = a.part("bracket", |m| {
    let b = m.profile(xz())
        .path(0.0, 0.0).line_to(44.0, 0.0).line_to(44.0, 5.0)
        .line_to(24.0, 5.0).line_to(0.0, 14.0).close()
        .sketch("wedge")
        .extrude(blind(24.0));

    // The tap: the inclined face has slope (24, -9) in xz, so its
    // inward normal is (-9, 0, -24)/|…| (exact, not rounded). The face
    // passes z = 9.5 at x = 12; 9 mm deep leaves 1 mm above the base.
    let n_len = sqrt(9.0 * 9.0 + 24.0 * 24.0);
    b = b.tapped_hole(
        axis_at([12.0, -12.0, 9.5], [-9.0 / n_len, 0.0, -24.0 / n_len]),
        "M6",
        blind(9.0),
        #{ clearance: 0.2, label: "angled_tap" });

    // The flange's mounting holes, countersunk in its top face.
    b = b.clearance_hole(
        axis_at([34.0, -5.5, 5.0], [0.0, 0.0, -1.0]),
        "M5",
        through_all(),
        #{ countersink: 0.2, label: "mount_a" });
    b = b.clearance_hole(
        axis_at([34.0, -18.5, 5.0], [0.0, 0.0, -1.0]),
        "M5",
        through_all(),
        #{ countersink: 0.2, label: "mount_b" });
});

// ISO 7089 M6 plain washer: 6.4 bore, 12 across, 1.6 thick.
let washer = a.part("washer", |m| {
    let bore = m.tool(|m| m.profile(xy_at(-0.5)).circle(3.2).sketch("bore").extrude(blind(2.6)));
    m.profile(xy()).circle(6.0).sketch("washer").extrude(blind(1.6)).subtract(bore);
});
let cap_screw = a.part("cap_screw", |m| {
    m.bolt("M6", 8.0, #{ head: "socket", thread_length: 8.0, label: "cap_screw" });
});
let flat_head = a.part("flat_head", |m| {
    m.bolt("M5", 10.0, #{ head: "countersunk", thread_length: 10.0, label: "flat_head" });
});

a.place(bracket);

// The tap's frame: z down the hole; x the world axis the tap seeds
// its thread phase from (the one most orthogonal to the axis, here
// y); y = z × x.
let n_len = sqrt(9.0 * 9.0 + 24.0 * 24.0);
let z = [-9.0 / n_len, 0.0, -24.0 / n_len];
let x = [0.0, 1.0, 0.0];
let y = [-z[2], 0.0, z[0]];
let entry = [12.0, -12.0, 9.5];

// The washer flat on the inclined face, its bore on the tap's axis.
a.place(washer).frame(entry, [-z[0], -z[1], -z[2]], x);

// The cap screw on the washer: lifted 1.6 mm up the axis, so turned
// back 1.6 turns (−576°, i.e. +144°) to stay in phase with the tap's
// thread (M6 lead 1). cos and sin of 144° are written exactly, so
// the native and browser builds place it bit-identically.
let lift = 1.6;
let c = -(1.0 + sqrt(5.0)) / 4.0;
let s = sqrt((5.0 - sqrt(5.0)) / 8.0);
a.place(cap_screw, "angled_screw").frame(
    [entry[0] - lift * z[0], entry[1] - lift * z[1], entry[2] - lift * z[2]],
    z,
    [c * x[0] + s * y[0], c * x[1] + s * y[1], c * x[2] + s * y[2]]);

// The flat-head screws seated flush in the countersinks, their tips
// standing 5 mm proud of the base where they would bite the bench.
a.place(flat_head, "mount_screw_a")
    .frame([34.0, -5.5, 5.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]);
a.place(flat_head, "mount_screw_b")
    .frame([34.0, -18.5, 5.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]);

a.assert_watertight();
a
