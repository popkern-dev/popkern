// A three-blade marine propeller. Each blade is LOFTED through four
// stadium sections this script computes — chord and pitch angle
// change along the span, so the blade twists like the real thing —
// then rotated from its build axis out to radial, phased 120° apart,
// and unioned onto a revolved hub with a tapered nose.
//
// NOT in the gallery yet: the blade ∪ hub union is the PK-443
// acceptance fixture (a chart-wrapping intersection loop on a closed
// lofted wall — the face split does not serve that class yet, and the
// union honestly declines). When PK-443 lands: run this script, add
// root fillets, and promote it to the gallery table
// (crates/popkern-cli/tests/gallery.rs).
let m = model();

let pi = 3.141592653589793;

// ---- Hub with nose cone: one revolved section ----
let hub = m.profile(xz())
    .path(0.0, -0.22).line_to(0.16, -0.22).line_to(0.16, 0.12)
    .line_to(0.03, 0.34).line_to(0.0, 0.34).close()
    .sketch("hub_section")
    .revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle(),
             #{ label: "hub" });

// ---- One blade: four computed sections, lofted, turned radial ----
// A section is a stadium (two lines, two semicircular edge caps —
// tangent-smooth all the way round, like a real blade's rounded
// leading and trailing edges) with chord c and thickness t, rotated
// in-plane by the local pitch angle. Sections stack along +z (the
// build axis); the finished blade is rotated about x so the span
// points +y, then phased about z.
fn blade_section(m, z, c, t, pitch_deg, label) {
    let pi = 3.141592653589793;
    let pitch = pitch_deg * pi / 180.0;
    let t2 = 0.5 * t;
    let c2 = 0.5 * c - t2;
    let rot = |x, y| [x * cos(pitch) - y * sin(pitch), x * sin(pitch) + y * cos(pitch)];
    let p0 = rot.call(-c2, -t2);
    let p1 = rot.call(c2, -t2);
    let v1 = rot.call(c2 + t2, 0.0);
    let p2 = rot.call(c2, t2);
    let p3 = rot.call(-c2, t2);
    let v2 = rot.call(-c2 - t2, 0.0);
    m.profile(xy_at(z))
        .path(p0[0], p0[1])
        .line_to(p1[0], p1[1])
        .arc_to(v1, p2)
        .line_to(p3[0], p3[1])
        .arc_to(v2, p0)
        .close()
        .sketch(label)
}

fn blade(m, phase_deg, tag) {
    let s0 = blade_section(m, 0.02, 0.15, 0.050, 52.0, tag + "_root");
    let s1 = blade_section(m, 0.13, 0.20, 0.034, 40.0, tag + "_mid1");
    let s2 = blade_section(m, 0.24, 0.17, 0.022, 30.0, tag + "_mid2");
    let s3 = blade_section(m, 0.34, 0.09, 0.012, 24.0, tag + "_tip");
    let b = m.loft([s0, s1, s2, s3], #{ label: tag + "_loft" });
    // Stand the blade out radially, then phase it about the shaft.
    let b = b.rotate(axis_at([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]), -90.0,
                     #{ label: tag + "_radial" });
    b.rotate(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), phase_deg,
             #{ label: tag + "_phase" })
}

let prop = hub;
let b1 = m.tool(|m| blade(m, 0.0, "blade_a"));
prop = prop.union(b1, #{ label: "join_a" });
let b2 = m.tool(|m| blade(m, 120.0, "blade_b"));
prop = prop.union(b2, #{ label: "join_b" });
let b3 = m.tool(|m| blade(m, 240.0, "blade_c"));
prop = prop.union(b3, #{ label: "join_c" });

m.assert_watertight();
m
