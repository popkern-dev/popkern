// A six-axis arm (UR-style kinematics, 850 mm reach) as an ASSEMBLY
// (design docs/architecture/six-axis-arm.md). Six cycloidal joint
// modules — the same modules as cycloidal_joint.rhai, sizes L, L, M,
// S, S, S — drive seven links: base, shoulder, upper arm, forearm,
// wrist 1, wrist 2 and an ISO 9409-1-50-4-M6 tool flange carrying a
// parallel-jaw gripper.
//
// Every link is authored at home (every joint at 0) in world
// coordinates and placed at identity, so its placement IS its motion
// from home and each joint axis is the same literal in both links'
// frames. A link is one rigid sub-assembly: its casting, the output
// side of the joint before it (flange, output pins, inner race) and
// the housing side of the joint after it (housing, stator, windings,
// ring pins, outer race), with their bolts. The rotor, the two cycloid
// discs and the roller cage of each joint are free roots that follow
// the joint's slider at their ratios.
//
// Sliders: j1..j6 (degrees, §4.1 limits), the tube lengths upper_arm
// and forearm, and grip (0..80). A tube length drives its tube part and
// the distance mate between the arm's two end castings, which are two
// roots joined along the tube's axis. The gripper's two fingers are one
// part, the second instance turned 180° about the tool axis; each
// slides on the gripper body, its jaw face grip/2 from the axis.
//
// Fasteners: socket cap screws (fastener family) in counterbored
// holes. Each mount ring clamps a housing's lip; each foot bolts to a
// joint's output flange. Every screw has a real thread and threads into
// a real tapped hole (PK-876): the housings' lips, the output flanges
// and the tool flange's four ISO M6 holes (a 0.1 fit clearance, 0.15
// for the M6: body_volume over-reads a 0.08-0.1 M6 flank, PK-901).
// Each tapped hole is one cut of a cached tool (~0.2-0.4 s on these
// castings). A screw mates a tap cut in its own insertion frame or an
// axial offset of whole leads (PK-476); `place_bolts` turns each screw
// by the leftover fraction of a lead between its head and the tapped
// face.
//
let a = assembly("six_axis_arm");
let pi = 3.141592653589793;

// ---- Size table (mm) ----
// n ring pins on radius big_r (pin radius rp), eccentricity e, backlash
// b; six output pins of radius op_r on radius rop; d = bearing roller.
let sizes = [
    #{ name: "L", od: 110.0, len: 118.0, pcd: 84.0, bolts: 12, tap: 4.2, screw: "M5",
       n: 30, big_r: 40.0, rp: 2.5, e: 0.9, b: 0.02, disc_t: 6.0, cam_r: 11.0,
       rop: 26.0, op_r: 3.0, shaft_r: 8.0, d: 5.0 },
    #{ name: "M", od: 92.0, len: 104.0, pcd: 70.0, bolts: 12, tap: 3.3, screw: "M4",
       n: 26, big_r: 33.0, rp: 2.0, e: 0.75, b: 0.02, disc_t: 5.0, cam_r: 9.0,
       rop: 21.5, op_r: 2.5, shaft_r: 6.5, d: 4.0 },
    #{ name: "S", od: 74.0, len: 88.0, pcd: 56.0, bolts: 8, tap: 3.3, screw: "M4",
       n: 22, big_r: 26.0, rp: 1.6, e: 0.6, b: 0.02, disc_t: 4.0, cam_r: 7.0,
       rop: 16.0, op_r: 1.8, shaft_r: 5.0, d: 3.0 },
];

// ---- Derived layout of one module (local frame: joint axis +Z) ----
fn layout(s) {
    let ro = s.od / 2.0;
    let ri = ro - 5.0;
    let len = s.len;
    let z_d0 = 0.5 * len;                       // disc A bottom
    let z_db = z_d0 + s.disc_t + 1.0;            // disc B bottom
    let z_s0 = z_d0 - 2.0;                       // ring-pin seat
    let z_s1 = z_db + s.disc_t + 2.0;
    let rim_in = ro - 10.0;
    let rb = rim_in - 2.0;                       // bearing pitch radius
    let half = (s.d + 0.2) / sqrt(2.0);          // roller diamond half-diagonal
    let hub_r = rb - half - 2.0;
    let mag_in = s.cam_r + 0.6 * (ri - s.cam_r) - 8.0;
    #{
        ro: ro, ri: ri, len: len, rim_in: rim_in,
        rs: s.big_r + s.rp + 0.01,               // seat bore: pins clear by 0.01
        z_d0: z_d0, z_db: z_db, z_s0: z_s0, z_s1: z_s1,
        floor: 4.0, dome: 0.18 * ro,
        rb: rb, zb: len - 16.5, half: half, hub_r: hub_r, hub_ri: hub_r - 2.0,
        flange_r: s.pcd / 2.0 + s.tap / 2.0 + 3.0,
        mag_in: mag_in, mag_out: mag_in + 5.0,
        st_in: mag_in + 5.8, st_out: ri - 0.5,
        z_st0: 12.0, z_st1: z_d0 - 8.0,
    }
}

// ---- Cycloid disc profile and its biarc fit ----
// c = #{ big_r, rr (pin radius + backlash), e, n }.
fn cyc_pt(c, t) {
    let k = 1.0 - c.n;
    let psi = atan(sin(k * t), c.big_r / (c.e * c.n) - cos(k * t));
    let x = c.big_r * cos(t) - c.rr * cos(t + psi) - c.e * cos(c.n * t);
    let y = -c.big_r * sin(t) + c.rr * sin(t + psi) + c.e * sin(c.n * t);
    [x, y]
}
fn cyc_tan(c, t) {
    let h = 1e-6;
    let p = cyc_pt(c, t - h);
    let q = cyc_pt(c, t + h);
    let dx = q[0] - p[0];
    let dy = q[1] - p[1];
    let l = sqrt(dx * dx + dy * dy);
    [dx / l, dy / l]
}
// The arc leaving p along tv and ending at q: [centre, radius], or []
// when it degenerates to the chord.
fn arc_centre(p, tv, q) {
    let nx = -tv[1];
    let ny = tv[0];
    let wx = q[0] - p[0];
    let wy = q[1] - p[1];
    let den = 2.0 * (wx * nx + wy * ny);
    let ww = wx * wx + wy * wy;
    if den.abs() < 1e-12 * sqrt(ww) { return []; }
    let s = ww / den;
    [[p[0] + s * nx, p[1] + s * ny], s.abs()]
}
// The biarc from (p0, t0) to (p1, t1): [joint, arc 1, arc 2].
fn biarc(p0, t0, p1, t1) {
    let vx = p1[0] - p0[0];
    let vy = p1[1] - p0[1];
    let vt = vx * (t0[0] + t1[0]) + vy * (t0[1] + t1[1]);
    let vv = vx * vx + vy * vy;
    let a2 = 2.0 * (t0[0] * t1[0] + t0[1] * t1[1] - 1.0);
    let d = 0.0;
    if a2.abs() < 1e-12 {
        d = vv / (4.0 * vt);
    } else {
        let disc = vt * vt - a2 * vv;
        d = (vt - sqrt(disc)) / a2;
        if d <= 0.0 { d = (vt + sqrt(disc)) / a2; }
    }
    let j = [0.5 * (p0[0] + d * t0[0] + p1[0] - d * t1[0]),
             0.5 * (p0[1] + d * t0[1] + p1[1] - d * t1[1])];
    [j, arc_centre(p0, t0, j), arc_centre(p1, t1, j)]
}
fn arc_dist(x, arc, p, q) {
    if arc.len() == 0 {
        let dx = q[0] - p[0];
        let dy = q[1] - p[1];
        let cr = (x[0] - p[0]) * dy - (x[1] - p[1]) * dx;
        return cr.abs() / sqrt(dx * dx + dy * dy);
    }
    let c = arc[0];
    let dx = x[0] - c[0];
    let dy = x[1] - c[1];
    (sqrt(dx * dx + dy * dy) - arc[1]).abs()
}
// The point of the minor arc between a and b (an arc_to `via`).
fn arc_mid(arc, a, b) {
    let mx = 0.5 * (a[0] + b[0]);
    let my = 0.5 * (a[1] + b[1]);
    if arc.len() == 0 { return [mx, my]; }
    let c = arc[0];
    let dx = mx - c[0];
    let dy = my - c[1];
    let l = sqrt(dx * dx + dy * dy);
    [c[0] + arc[1] * dx / l, c[1] + arc[1] * dy / l]
}
// Biarcs over [t0, t1], halved until every sample is within tol:
// [segments ([via, to] pairs), worst sampled deviation].
fn fit_span(c, t0, t1, tol, depth) {
    let p0 = cyc_pt(c, t0);
    let p1 = cyc_pt(c, t1);
    let bi = biarc(p0, cyc_tan(c, t0), p1, cyc_tan(c, t1));
    let j = bi[0];
    let worst = 0.0;
    for k in 1..16 {
        let x = cyc_pt(c, t0 + (t1 - t0) * k / 16.0);
        let d1 = arc_dist(x, bi[1], p0, j);
        let d2 = arc_dist(x, bi[2], j, p1);
        let d = if d1 < d2 { d1 } else { d2 };
        if d > worst { worst = d; }
    }
    if worst > tol && depth < 12 {
        let tm = 0.5 * (t0 + t1);
        let l = fit_span(c, t0, tm, tol, depth + 1);
        let r = fit_span(c, tm, t1, tol, depth + 1);
        let segs = l[0];
        segs += r[0];
        return [segs, if l[1] > r[1] { l[1] } else { r[1] }];
    }
    [[[arc_mid(bi[1], p0, j), j], [arc_mid(bi[2], j, p1), p1]], worst]
}
// The whole flank: [segments, worst deviation]. Four spans per lobe.
fn cycloid_fit(s, pi) {
    let c = #{ big_r: s.big_r, rr: s.rp + s.b, e: s.e, n: s.n };
    let spans = 4 * (s.n - 1);
    let segs = [];
    let worst = 0.0;
    for k in 0..spans {
        let f = fit_span(c, 2.0 * pi * k / spans, 2.0 * pi * (k + 1) / spans, s.b / 4.0, 0);
        segs += f[0];
        if f[1] > worst { worst = f[1]; }
    }
    [segs, worst, cyc_pt(c, 0.0)]
}

// ---- Part bodies (each a function: closures stay shallow) ----
// Size s's bolt circle of radius ring_r: real tapped holes (PK-876)
// from the face z = entry down -z, `depth` deep, with a 0.1 fit
// clearance. Every hole phases its thread from the local +x seed; the
// screws are placed to match (`place_bolts`).
fn bolt_taps(body, s, ring_r, entry, depth, pi) {
    for k in 0..s.bolts {
        let ang = 2.0 * pi * k / s.bolts;
        body = body.tapped_hole(axis_at([ring_r * cos(ang), ring_r * sin(ang), entry], [0.0, 0.0, -1.0]),
            s.screw, blind(depth), #{ clearance: 0.1, label: "tapped_" + s.screw + "_" + k });
    }
    body
}
// The housing's outer solid: a revolved cylinder with a domed cap and
// a mounting lip at the output end (the bolt circle a casting clamps).
fn housing_solid(m, g) {
    let ro = g.ro;
    let zc = (ro * ro - g.dome * g.dome) / (2.0 * g.dome);
    let rho = zc + g.dome;
    let vr = 0.7 * ro;
    let vz = zc - sqrt(rho * rho - vr * vr);
    let p = m.profile(xz()).path(0.0, -g.dome).arc_to([vr, vz], [ro, 0.0]);
    p = p.line_to(ro, g.len - 8.0).line_to(ro + 9.0, g.len - 8.0).line_to(ro + 9.0, g.len);
    p = p.line_to(0.0, g.len).close();
    p.sketch("shell").revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle())
}
// The cavity, revolved: floor, the ring-pin seat, the bolt rim.
fn housing_cavity(t, g) {
    let p = t.profile(xz()).path(0.0, g.floor).line_to(g.ri, g.floor).line_to(g.ri, g.z_s0);
    p = p.line_to(g.rs, g.z_s0).line_to(g.rs, g.z_s1).line_to(g.ri, g.z_s1);
    p = p.line_to(g.ri, g.len - 10.0).line_to(g.rim_in, g.len - 10.0).line_to(g.rim_in, g.len + 1.0);
    p = p.line_to(0.0, g.len + 1.0).close();
    p.sketch("cavity").revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle())
}
// 24 radial fins as ONE polygon whose inner chords lie inside the cup.
fn fin_ring(t, g, pi) {
    let pts = [];
    let w = 1.5;
    let root = g.ro - 2.0;
    let tip = g.ro + 8.0;
    for k in 0..24 {
        let ang = 2.0 * pi * k / 24.0;
        let ca = cos(ang);
        let sa = sin(ang);
        pts.push([root * ca + w * sa, root * sa - w * ca]);
        pts.push([tip * ca + w * sa, tip * sa - w * ca]);
        pts.push([tip * ca - w * sa, tip * sa + w * ca]);
        pts.push([root * ca - w * sa, root * sa + w * ca]);
    }
    t.profile(xy_at(0.2 * g.len)).polygon(pts).sketch("fins").extrude(blind(0.65 * g.len))
}
// A flat pad tangent-plus-1 mm on the wall below the fins, opposite
// the cable boss: a rectangle on a plane through the wall, extruded
// outward.
fn label_pad(t, g, pi) {
    let ang = pi + pi / 24.0;
    let n = [cos(ang), sin(ang), 0.0];
    let depth = g.ro - 2.0;
    let origin = [depth * n[0], depth * n[1], 0.1 * g.len];
    t.profile(plane_at(origin, n)).rect(0.08 * g.len, 0.35 * g.ro).sketch("pad")
        .extrude(blind(3.0))
}
fn radial_rod(t, z, ang, r, length, label) {
    t.profile(plane_at([0.0, 0.0, z], [cos(ang), sin(ang), 0.0])).circle(r).sketch(label)
        .extrude(blind(length))
}
// The fins are ONE polygon (a filled star): unioned onto the SOLID
// body, its centre is absorbed; the cavity is cut after.
fn housing_body(m, g, s, pi) {
    let body = housing_solid(m, g);
    body = body.union(m.tool(|t| fin_ring(t, g, pi)), #{ label: "fins" });
    body = body.fillet(edges().straight().direction("z").between("plane", "cylinder"), 1.0,
        #{ label: "fin_roots" });
    // Cable boss between the dome and the fins.
    let zb = 0.1 * g.len;
    let ang = pi / 24.0;
    body = body.union(m.tool(|t| radial_rod(t, zb, ang, 5.0, g.ro + 6.0, "cable_boss")),
        #{ label: "cable_boss" });
    // Label pad opposite the boss: a flat raised 1 mm off the wall.
    body = body.union(m.tool(|t| label_pad(t, g, pi)), #{ label: "label_pad" });
    body = body.subtract(m.tool(|t| housing_cavity(t, g)), #{ label: "cavity" });
    body = body.subtract(m.tool(|t| radial_rod(t, zb, ang, 2.5, g.ro + 8.0, "cable_bore")),
        #{ label: "cable_bore" });
    // Bolt circle through the 8 mm mounting lip: real tapped holes
    // (PK-876) whose thread runs out of the lip's underside — 10.5
    // deep, the floor in the air between the lip and the fins (7.6
    // below it). Ending 9 deep, the M4 thread's tip turned 0.12 past
    // the underside and the cut declined.
    bolt_taps(body, s, g.ro + 4.5, g.len, 10.5, pi)
}
// The slotted bore of the stator: bore arcs and 12 rectangular slots.
fn stator_star(g, slot_w, slot_r, pi) {
    let star = [];
    let hw = slot_w / 2.0;
    let a_in = asin(hw / g.st_in);
    for k in 0..12 {
        let a0 = 2.0 * pi * k / 12.0;
        let ca = cos(a0);
        let sa = sin(a0);
        let a_prev = a0 - 2.0 * pi / 12.0 + a_in;
        for q in 1..4 {
            let t = a_prev + (a0 - a_in - a_prev) * q / 4.0;
            star.push([g.st_in * cos(t), g.st_in * sin(t)]);
        }
        star.push([g.st_in * cos(a0 - a_in), g.st_in * sin(a0 - a_in)]);
        star.push([slot_r * ca + hw * sa, slot_r * sa - hw * ca]);
        star.push([slot_r * ca - hw * sa, slot_r * sa + hw * ca]);
        star.push([g.st_in * cos(a0 + a_in), g.st_in * sin(a0 + a_in)]);
    }
    star
}
fn stator_body(m, g, slot_w, slot_r, pi) {
    let ring = m.profile(xy_at(g.z_st0)).circle(g.st_out).sketch("stack")
        .extrude(blind(g.z_st1 - g.z_st0));
    let star = stator_star(g, slot_w, slot_r, pi);
    let slots = m.tool(|t| t.profile(xy_at(g.z_st0 - 1.0)).polygon(star).sketch("slots")
        .extrude(blind(g.z_st1 - g.z_st0 + 2.0)));
    ring.subtract(slots, #{ label: "slots" })
}
fn winding_body(m, g, slot_w, slot_r) {
    let r0 = g.st_in + 0.8;
    m.profile(xy_at(g.z_st0 - 1.5)).rect(slot_r - r0 - 0.3, slot_w - 0.6)
        .at(0.5 * (slot_r + r0 - 0.3), 0.0).sketch("coil").extrude(blind(g.z_st1 - g.z_st0 + 3.0))
}
fn rod(t, z, cx, cy, r, length, label) {
    t.profile(xy_at(z)).circle(r).at(cx, cy).sketch(label).extrude(blind(length))
}
// A shaft with the magnet hub and two cams 180° apart.
fn shaft_body(m, g, s) {
    let z0 = g.floor + 0.5;
    let core = rod(m, z0, 0.0, 0.0, s.shaft_r, g.z_s1 + 2.0 - z0, "shaft");
    core = core.union(m.tool(|t| rod(t, g.z_st0, 0.0, 0.0, g.mag_in, g.z_st1 - g.z_st0, "magnet_hub")),
        #{ label: "magnet_hub" });
    core = core.union(m.tool(|t| rod(t, g.z_d0, s.e, 0.0, s.cam_r, s.disc_t, "cam_a")),
        #{ label: "cam_a" });
    core.union(m.tool(|t| rod(t, g.z_db, -s.e, 0.0, s.cam_r, s.disc_t, "cam_b")),
        #{ label: "cam_b" })
}
fn magnet_body(m, g, pi) {
    let h = pi / 10.0 - 0.05;
    let r0 = g.mag_in;
    let r1 = g.mag_out;
    let p = m.profile(xy_at(g.z_st0 + 1.0)).path(r0 * cos(-h), r0 * sin(-h));
    p = p.line_to(r1 * cos(-h), r1 * sin(-h)).arc_to([r1, 0.0], [r1 * cos(h), r1 * sin(h)]);
    p = p.line_to(r0 * cos(h), r0 * sin(h)).arc_to([r0, 0.0], [r0 * cos(-h), r0 * sin(-h)]).close();
    p.sketch("segment").extrude(blind(g.z_st1 - g.z_st0 - 2.0))
}
// A 64-window wheel, the windows open at the rim: one sketch.
fn encoder_body(m, g, s, pi) {
    let r1 = g.mag_in;
    let r0 = r1 - 3.0;
    let hw = 0.25 * 2.0 * pi / 64.0;
    let pts = [];
    for k in 0..64 {
        let t0 = 2.0 * pi * k / 64.0;
        pts.push([r1 * cos(t0 - hw), r1 * sin(t0 - hw)]);
        pts.push([r0 * cos(t0 - hw), r0 * sin(t0 - hw)]);
        pts.push([r0 * cos(t0 + hw), r0 * sin(t0 + hw)]);
        pts.push([r1 * cos(t0 + hw), r1 * sin(t0 + hw)]);
    }
    let wheel = m.profile(xy_at(g.floor + 1.5)).polygon(pts).sketch("wheel").extrude(blind(1.5));
    wheel.subtract(m.tool(|t| rod(t, g.floor, 0.0, 0.0, s.shaft_r + 0.05, 4.0, "hub_bore")),
        #{ label: "hub_bore" })
}
// The cam bore and six output-pin windows (pin radius + e + 0.05).
fn disc_cut(t, s, pi) {
    let all = rod(t, -1.0, 0.0, 0.0, s.cam_r + 0.05, s.disc_t + 2.0, "cam_bore");
    let window = s.op_r + s.e + 0.05;
    for k in 0..6 {
        let ang = 2.0 * pi * k / 6.0;
        all = all.union(rod(t, -1.0, s.rop * cos(ang), s.rop * sin(ang), window, s.disc_t + 2.0,
            "window_" + k));
    }
    all
}
fn disc_body(m, s, fit, pi) {
    let p = m.profile(xy_at(0.0)).path(fit[2][0], fit[2][1]);
    for seg in fit[0] { p = p.arc_to(seg[0], seg[1]); }
    let blank = p.sketch("flank").extrude(blind(s.disc_t));
    blank.subtract(m.tool(|t| disc_cut(t, s, pi)), #{ label: "bores" })
}
// Cross-roller races around a diamond groove (rollers at ±45°).
fn race_body(m, g, outer) {
    let h = g.half + 1.5;
    let gap = 0.1;
    let p = m.profile(xz());
    if outer {
        p = p.path(g.rb + gap, g.zb + h).line_to(g.ri - 0.1, g.zb + h).line_to(g.ri - 0.1, g.zb - h);
        p = p.line_to(g.rb + gap, g.zb - h).line_to(g.rb + gap, g.zb - g.half + gap);
        p = p.line_to(g.rb + g.half, g.zb).line_to(g.rb + gap, g.zb + g.half - gap).close();
    } else {
        p = p.path(g.hub_r, g.zb - h).line_to(g.rb - gap, g.zb - h);
        p = p.line_to(g.rb - gap, g.zb - g.half + gap).line_to(g.rb - g.half, g.zb);
        p = p.line_to(g.rb - gap, g.zb + g.half - gap).line_to(g.rb - gap, g.zb + h);
        p = p.line_to(g.hub_r, g.zb + h).close();
    }
    p.sketch("race").revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle())
}
// The output flange: a hub down to the inner race and a top plate
// with the output bolt circle (real tapped holes, PK-876).
fn flange_body(m, g, s, pi) {
    let top = g.len + 1.0;
    let z0 = g.zb - g.half - 1.5;
    let p = m.profile(xz()).path(g.hub_ri, z0).line_to(g.hub_r, z0).line_to(g.hub_r, top);
    p = p.line_to(g.flange_r, top).line_to(g.flange_r, top + 10.0).line_to(0.0, top + 10.0);
    p = p.line_to(0.0, top).line_to(g.hub_ri, top).close();
    let body = p.sketch("flange").revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle());
    bolt_taps(body, s, s.pcd / 2.0, top + 10.0, 8.0, pi)
}

// ---- One module's part definitions ----
fn module_parts(a, s, pi) {
    let g = layout(s);
    let fit = cycloid_fit(s, pi);
    if fit[1] > s.b / 4.0 {
        throw "disc " + s.name + ": backlash " + s.b + " is below the profile fit deviation "
            + fit[1] * 4.0 + " (the fit must stay within b/4)";
    }
    print("disc " + s.name + ": " + fit[0].len() + " arcs, fit deviation " + fit[1]
        + " (b/4 = " + s.b / 4.0 + ")");
    let slot_w = 0.45 * 2.0 * pi * g.st_in / 12.0;
    let slot_r = g.st_in + 0.5 * (g.st_out - g.st_in);
    let n = s.name;
    let p = #{ g: g };
    p.housing = a.part("joint_housing_" + n, |m| housing_body(m, g, s, pi));
    p.stator = a.part("stator_" + n, |m| stator_body(m, g, slot_w, slot_r, pi));
    p.winding = a.part("winding_" + n, |m| winding_body(m, g, slot_w, slot_r));
    p.shaft = a.part("shaft_" + n, |m| shaft_body(m, g, s));
    p.magnet = a.part("magnet_" + n, |m| magnet_body(m, g, pi));
    p.encoder = a.part("encoder_" + n, |m| encoder_body(m, g, s, pi));
    p.disc = a.part("disc_" + n, |m| disc_body(m, s, fit, pi));
    p.ring_pin = a.part("ring_pin_" + n, |m| rod(m, g.z_s0 + 0.5, s.big_r, 0.0, s.rp,
        g.z_s1 - g.z_s0 - 1.0, "pin"));
    p.outer_race = a.part("bearing_outer_" + n, |m| race_body(m, g, true));
    p.inner_race = a.part("bearing_inner_" + n, |m| race_body(m, g, false));
    p.roller = a.part("roller_" + n, |m| rod(m, -s.d / 2.0, 0.0, 0.0, s.d / 2.0, s.d, "roller"));
    p.flange = a.part("output_" + n, |m| flange_body(m, g, s, pi));
    p.output_pin = a.part("output_pin_" + n, |m| rod(m, g.z_d0 - 1.0, s.rop, 0.0, s.op_r,
        g.len + 1.0 - (g.z_d0 - 1.0), "pin"));
    let looks = [
        ["housing", "aluminium", 0.7, 0.35], ["stator", "gray", 0.5, 0.5],
        ["winding", "copper", 0.8, 0.3], ["shaft", "steel", 0.85, 0.25],
        ["magnet", "red", 0.4, 0.4], ["encoder", "black", 0.2, 0.6],
        ["disc", "steel", 0.85, 0.2], ["ring_pin", "steel", 0.9, 0.15],
        ["outer_race", "steel", 0.9, 0.2], ["inner_race", "steel", 0.9, 0.2],
        ["roller", "steel", 0.9, 0.2], ["flange", "blue", 0.6, 0.35],
        ["output_pin", "steel", 0.9, 0.15],
    ];
    for l in looks {
        a.appearance(p[l[0]], #{ color: l[1], metal: l[2], rough: l[3] });
    }
    p
}

// ---- Rigid transforms (4x4, row-major) ----
fn m_ident() { [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]] }
fn m_mul(a, b) {
    let r = m_ident();
    for i in 0..4 {
        for j in 0..4 {
            let s = 0.0;
            for k in 0..4 { s += a[i][k] * b[k][j]; }
            r[i][j] = s;
        }
    }
    r
}
fn m_rz(t) { let c = cos(t); let s = sin(t); [[c, -s, 0.0, 0.0], [s, c, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]] }
fn m_rx(t) { let c = cos(t); let s = sin(t); [[1.0, 0.0, 0.0, 0.0], [0.0, c, -s, 0.0], [0.0, s, c, 0.0], [0.0, 0.0, 0.0, 1.0]] }
fn m_tr(x, y, z) { [[1.0, 0.0, 0.0, x], [0.0, 1.0, 0.0, y], [0.0, 0.0, 1.0, z], [0.0, 0.0, 0.0, 1.0]] }
fn m_inv(m) {
    let r = m_ident();
    for i in 0..3 { for j in 0..3 { r[i][j] = m[j][i]; } }
    for i in 0..3 { r[i][3] = -(r[i][0] * m[0][3] + r[i][1] * m[1][3] + r[i][2] * m[2][3]); }
    r
}
fn m_point(m, p) {
    [m[0][0] * p[0] + m[0][1] * p[1] + m[0][2] * p[2] + m[0][3],
     m[1][0] * p[0] + m[1][1] * p[1] + m[1][2] * p[2] + m[1][3],
     m[2][0] * p[0] + m[2][1] * p[1] + m[2][2] * p[2] + m[2][3]]
}
fn m_dir(m, d) {
    [m[0][0] * d[0] + m[0][1] * d[1] + m[0][2] * d[2],
     m[1][0] * d[0] + m[1][1] * d[1] + m[1][2] * d[2],
     m[2][0] * d[0] + m[2][1] * d[1] + m[2][2] * d[2]]
}
// One standard Denavit-Hartenberg step.
fn m_dh(theta, d, a, alpha) { m_mul(m_mul(m_rz(theta), m_tr(0.0, 0.0, d)), m_mul(m_tr(a, 0.0, 0.0), m_rx(alpha))) }

// ---- Arm layout (§4.1; mm) ----
// Each joint's module: its size, and the offset c of its module origin
// (dome end) along the joint axis from the DH frame before it. Every
// module's output points along its DH axis z.
fn arm_layout(pi) {
    #{
        d: [162.5, 0.0, 0.0, 133.3, 99.7, 99.6],
        a: [0.0, -425.0, -392.0, 0.0, 0.0, 0.0],
        alpha: [pi / 2.0, 0.0, 0.0, pi / 2.0, -pi / 2.0, 0.0],
        limits: [[-360.0, 360.0], [-360.0, 360.0], [-170.0, 170.0],
                 [-360.0, 360.0], [-360.0, 360.0], [-360.0, 360.0]],
        sizes: [0, 0, 1, 2, 2, 2],
        c: [-50.0, -128.0, -102.0, -72.0, -58.5, -17.4],
        // Tube axes (y, at the joint height z = d1) and end-casting
        // reach along the arm, from the joint axes.
        upper_y: -45.0, upper_r: 42.5, upper_ri: 38.5, upper_root: 90.0, upper_end: 90.0,
        fore_y: -50.0, fore_r: 33.0, fore_ri: 29.0, fore_root: 80.0, fore_end: 77.0,
    }
}
// World transforms of the DH frames at home (every joint at 0):
// frames[k] for k = 0..6, frame 0 the base.
fn dh_frames(arm) {
    let frames = [m_ident()];
    for k in 0..6 {
        frames.push(m_mul(frames[k], m_dh(0.0, arm.d[k], arm.a[k], arm.alpha[k])));
    }
    frames
}
// Per joint: its module's world frame at home, the joint axis, and the
// casting dimensions its bolt circles set. A mount ring (on the link
// before the joint) clamps the housing lip; a foot (on the link after)
// bolts to the output flange.
fn joint_mods(arm, frames, sizes) {
    let mods = [];
    for k in 0..6 {
        let s = sizes[arm.sizes[k]];
        let g = layout(s);
        let big = s.name == "L";
        let cb_r = if big { 5.0 } else { 4.25 };
        let f = frames[k];
        mods.push(#{
            k: k, s: s, g: g,
            M: m_mul(f, m_tr(0.0, 0.0, arm.c[k])),
            o: [f[0][3], f[1][3], f[2][3]], z: [f[0][2], f[1][2], f[2][2]],
            ring_in: g.flange_r + 1.4, ring_out: g.ro + 12.0, rl: g.ro + 4.5,
            foot_r: s.pcd / 2.0 + cb_r + 3.0,
            clear_r: if big { 2.75 } else { 2.25 }, cb_r: cb_r, cb_d: if big { 6.0 } else { 5.0 },
            screw: if big { "M5" } else { "M4" }, lead: if big { 0.8 } else { 0.7 },
        });
    }
    mods
}

// ---- Placement by a rigid transform ----
fn place_at(inst, m) {
    inst.frame([m[0][3], m[1][3], m[2][3]], [m[0][2], m[1][2], m[2][2]], [m[0][0], m[1][0], m[2][0]])
}

// ---- Primitive solids in world coordinates (home pose) ----
fn rod_at(t, p, n, r, length, label) {
    t.profile(plane_at(p, n)).circle(r).sketch(label).extrude(blind(length))
}
fn box(t, x0, x1, y0, y1, z0, z1, label) {
    t.profile(xy_at(z0)).rect(x1 - x0, y1 - y0).at(0.5 * (x0 + x1), 0.5 * (y0 + y1))
        .sketch(label).extrude(blind(z1 - z0))
}
// A disc on module md's axis from ζ = z0 to z1.
fn axial_disc(t, md, z0, z1, r, label) {
    rod_at(t, m_point(md.M, [0.0, 0.0, z0]), m_dir(md.M, [0.0, 0.0, 1.0]), r, z1 - z0, label)
}
// Counterbored clearance holes on circle rc of module md: the face at
// ζ = zf, the pocket `depth` deep, the clearance hole to ζ = zf −
// through. Each pocket also opens `ext` outward past the face: the
// bolt's access bore through whatever casting lies behind the face.
fn cbore_tool(t, md, rc, count, zf, depth, through, ext, pi) {
    let all = ();
    let d = m_dir(md.M, [0.0, 0.0, -1.0]);
    for k in 0..count {
        let ang = 2.0 * pi * k / count;
        let u = rc * cos(ang);
        let v = rc * sin(ang);
        let pocket = rod_at(t, m_point(md.M, [u, v, zf + ext]), d, md.cb_r, ext + depth, "pocket_" + k);
        let hole = rod_at(t, m_point(md.M, [u, v, zf - depth + 0.5]), d, md.clear_r,
            through - depth + 0.5, "clear_" + k);
        let both = pocket.union(hole);
        all = if k == 0 { both } else { all.union(both) };
    }
    all
}
// The bolts of a cbore_tool circle, heads on the pocket floors, each
// threading into the tapped hole whose entry face is `gap` below its
// head. Thread mating phase (PK-476): a screw mates a tap cut in its
// own insertion frame, or one offset along the axis by whole leads; the
// head sits `gap` out of the hole, so the screw turns by the leftover
// fraction of a lead about the outward axis from the hole's +x seed.
fn place_bolts(sub, bolt, md, rc, count, zf, depth, gap, name, pi) {
    let d = m_dir(md.M, [0.0, 0.0, -1.0]);
    let leads = gap / md.lead;
    let phase = 2.0 * pi * (leads - leads.floor());
    let x = m_dir(md.M, [cos(phase), sin(phase), 0.0]);
    for k in 0..count {
        let ang = 2.0 * pi * k / count;
        let at = m_point(md.M, [rc * cos(ang), rc * sin(ang), zf - depth]);
        sub.place(bolt, name + "_" + k).frame(at, d, x);
    }
}
// The mount ring's and the foot's bolt circles (module md); `ext` is
// the access bore: past any casting behind the face, or short of it.
fn ring_holes(t, md, ext, pi) {
    cbore_tool(t, md, md.rl, md.s.bolts, md.g.len + 10.0, md.cb_d, 10.5, ext, pi)
}
fn foot_holes(t, md, ext, pi) {
    cbore_tool(t, md, md.s.pcd / 2.0, md.s.bolts, md.g.len + 23.0, md.cb_d, 12.5, ext, pi)
}
fn foot_disc(t, md) { axial_disc(t, md, md.g.len + 11.0, md.g.len + 23.0, md.foot_r, "foot") }
fn ring_bore(t, md) { axial_disc(t, md, md.g.len - 1.0, md.g.len + 11.0, md.ring_in, "bore") }
// ---- Link castings (world coordinates, home pose) ----
// Base: a mounting plate with four counterbored M8 slots, a column
// round the j1 housing ending in its mount ring, a connector boss.
// A stadium slot (centres 12 apart along x, width w) from z0 to z1,
// sketched as one profile.
fn slot_tool(t, x, y, w, z0, z1, label) {
    let h = w / 2.0;
    let p = t.profile(xy_at(z0)).path(x - 6.0, y - h).line_to(x + 6.0, y - h);
    p = p.arc_to([x + 6.0 + h, y], [x + 6.0, y + h]).line_to(x - 6.0, y + h);
    p = p.arc_to([x - 6.0 - h, y], [x - 6.0, y - h]).close();
    p.sketch(label).extrude(blind(z1 - z0))
}
fn base_slots(t) {
    let all = ();
    let k = 0;
    for c in [[95.0, 95.0], [-95.0, 95.0], [-95.0, -95.0], [95.0, -95.0]] {
        let one = slot_tool(t, c[0], c[1], 9.0, -81.0, -65.0, "slot_" + k)
            .union(slot_tool(t, c[0], c[1], 15.0, -74.5, -50.0, "slot_cb_" + k));
        all = if k == 0 { one } else { all.union(one) };
        k += 1;
    }
    all
}
fn base_body(m, m1, pi) {
    let body = m.profile(xy_at(-80.0)).rect(240.0, 240.0).sketch("plate").extrude(blind(14.0));
    let top = m1.M[2][3] + m1.g.len;
    let col = m.tool(|t| {
        let p = t.profile(xz()).path(66.0, -70.0).line_to(76.0, -70.0).line_to(76.0, top + 10.0);
        p = p.line_to(m1.ring_in, top + 10.0).line_to(m1.ring_in, top).line_to(66.0, top).close();
        p.sketch("column").revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle())
    });
    body = body.union(col, #{ label: "column" });
    body = body.union(m.tool(|t| rod_at(t, [70.0, 0.0, -30.0], [1.0, 0.0, 0.0], 14.0, 30.0, "boss")),
        #{ label: "connector_boss" });
    body = body.subtract(m.tool(|t| rod_at(t, [60.0, 0.0, -30.0], [1.0, 0.0, 0.0], 8.0, 41.0, "port")),
        #{ label: "connector_port" });
    body = body.subtract(m.tool(|t| base_slots(t)), #{ label: "slots_M8" });
    body.subtract(m.tool(|t| ring_holes(t, m1, 30.0, pi)), #{ label: "j1_mount" })
}
// Shoulder: a foot on j1's flange and an upright plate carrying j2's
// mount ring.
fn shoulder_body(m, m1, m2, pi) {
    let body = foot_disc(m, m1);
    body = body.union(m.tool(|t| axial_disc(t, m2, m2.g.len, m2.g.len + 10.0, m2.ring_out, "seat")),
        #{ label: "seat" });
    // The web sits inside the seat's band: no coplanar faces to merge.
    body = body.union(m.tool(|t| {
        let y0 = m_point(m2.M, [0.0, 0.0, m2.g.len + 10.0])[1];
        box(t, -30.0, 30.0, y0 + 1.0, y0 + 9.0, m1.M[2][3] + m1.g.len + 17.0, 135.0, "web")
    }), #{ label: "web" });
    body = body.subtract(m.tool(|t| ring_bore(t, m2)), #{ label: "j2_bore" });
    // Pockets stop short of the seat above them.
    body = body.subtract(m.tool(|t| foot_holes(t, m1, 15.0, pi)), #{ label: "j1_foot" });
    body.subtract(m.tool(|t| ring_holes(t, m2, 30.0, pi)), #{ label: "j2_mount" })
}
// An end casting on a joint's flange: a foot wide enough to carry a
// boss along the arm (x from x0 to x1) clear of its bolt pockets, which
// stay open to the flange face.
fn arm_root_body(m, md, ty, r, rf, x0, x1, pi) {
    let body = axial_disc(m, md, md.g.len + 11.0, md.g.len + 23.0, rf, "foot");
    let zc = md.o[2];
    body = body.union(m.tool(|t| rod_at(t, [x0, ty, zc], [1.0, 0.0, 0.0], r, x1 - x0, "boss")),
        #{ label: "boss" });
    body.subtract(m.tool(|t| foot_holes(t, md, 5.0, pi)), #{ label: "foot_bolts" })
}
// The far end casting: the next joint's mount ring, a web out to the
// tube (x to x1, half-height hw, inset in the ring's band: no coplanar
// faces to merge; its top 1 mm below the boss top, never tangent), and
// the boss the tube butts against (x from x0 to x1), starting 2 mm
// clear of the ring so its cap never grazes the rim.
// The ring and web are two solids, the web inset in the ring's band.
fn arm_end_body(m, md, ty, r, x0, x1, hw, pi) {
    let body = axial_disc(m, md, md.g.len, md.g.len + 10.0, md.ring_out, "ring");
    let c = m_point(md.M, [0.0, 0.0, md.g.len]);
    let y0 = c[1].min(m_point(md.M, [0.0, 0.0, md.g.len + 10.0])[1]);
    let xs = if x1 > c[0] { c[0] + md.ring_in } else { c[0] - md.ring_in };
    body = body.union(m.tool(|t| box(t, xs.min(x1), xs.max(x1), y0 + 1.0, y0 + 8.0, c[2] - hw, c[2] + hw,
        "web_sketch")), #{ label: "web" });
    body = body.union(m.tool(|t| rod_at(t, [x0, ty, md.o[2]], [1.0, 0.0, 0.0], r, x1 - x0, "boss")),
        #{ label: "boss" });
    body = body.subtract(m.tool(|t| ring_bore(t, md)), #{ label: "bore" });
    body.subtract(m.tool(|t| ring_holes(t, md, 5.0, pi)), #{ label: "mount_bolts" })
}
// A tube on the line (y, z) from x1 toward −x for `length` (a part
// parameter expression), outer r, inner ri: the end at x1 is fixed,
// the far end follows the parameter.
fn tube_body(m, x1, y, z, r, ri, length, bore_length) {
    let body = m.profile(plane_at([x1, y, z], [-1.0, 0.0, 0.0])).circle(r).sketch("tube")
        .extrude(blind(expr(length)));
    body.subtract(m.tool(|t| t.profile(plane_at([x1 + 1.0, y, z], [-1.0, 0.0, 0.0])).circle(ri)
        .sketch("bore").extrude(blind(expr(bore_length)))), #{ label: "bore" })
}
// Wrist 1: a foot on j4's flange, j5's mount ring, and a two-step web
// kept within 42 mm of the j4 axis (it turns past the forearm's boss).
fn wrist_1_body(m, m4, m5, pi) {
    let x = m4.o[0];
    let body = foot_disc(m, m4);
    body = body.union(m.tool(|t| axial_disc(t, m5, m5.g.len, m5.g.len + 10.0, m5.ring_out, "ring_sketch")),
        #{ label: "ring" });
    body = body.union(m.tool(|t| box(t, x - 14.0, x + 14.0, -86.0, -30.0, 126.0, 150.0, "web_sketch")),
        #{ label: "web" });
    body = body.union(m.tool(|t| box(t, x - 12.0, x + 12.0, -92.0, -40.0, 124.0, 129.0, "web_low_sketch")),
        #{ label: "web_low" });
    body = body.subtract(m.tool(|t| ring_bore(t, m5)), #{ label: "j5_bore" });
    body = body.subtract(m.tool(|t| foot_holes(t, m4, 55.0, pi)), #{ label: "j4_foot" });
    body.subtract(m.tool(|t| ring_holes(t, m5, 30.0, pi)), #{ label: "j5_mount" })
}
// Wrist 2: a foot under j5's flange, j6's mount ring, and a web that
// steps over j6's housing lip.
fn wrist_2_body(m, m5, m6, pi) {
    let x = m5.o[0];
    let body = foot_disc(m, m5);
    body = body.union(m.tool(|t| axial_disc(t, m6, m6.g.len, m6.g.len + 10.0, 54.0, "ring_sketch")),
        #{ label: "ring" });
    body = body.union(m.tool(|t| box(t, x - 18.0, x + 18.0, -205.5, -150.0, 110.5, 121.0, "web_sketch")),
        #{ label: "web" });
    body = body.union(m.tool(|t| box(t, x - 16.0, x + 16.0, -213.0, -204.5, 104.0, 121.5, "web_low_sketch")),
        #{ label: "web_low" });
    body = body.subtract(m.tool(|t| ring_bore(t, m6)), #{ label: "j6_bore" });
    body = body.subtract(m.tool(|t| foot_holes(t, m5, 10.0, pi)), #{ label: "j5_foot" });
    body.subtract(m.tool(|t| ring_holes(t, m6, 30.0, pi)), #{ label: "j6_mount" })
}
// Tool flange, ISO 9409-1-50-4-M6: 4 REAL M6 tapped holes (the tool
// interface, PK-801) and a Ø6 dowel hole on PCD 50, a Ø31.5 centring
// bore; its own bolts sit in counterbores from the flange face.
fn iso_taps(body, md, face, pi) {
    let d = m_dir(md.M, [0.0, 0.0, -1.0]);
    for k in 0..4 {
        let ang = pi / 8.0 + pi * k / 2.0;
        body = body.tapped_hole(axis_at(m_point(md.M, [25.0 * cos(ang), 25.0 * sin(ang), face]), d),
            "M6", blind(10.0), #{ clearance: 0.15, label: "tapped_M6_" + k });
    }
    body
}
// A hole's insertion-frame x seed along d: the world axis most
// orthogonal to it (popkern_fastener::holes::insertion_x_seed).
fn insertion_seed(d) {
    let x = abs(d[0]);
    let y = abs(d[1]);
    let z = abs(d[2]);
    if x <= y && x <= z { [1.0, 0.0, 0.0] } else if y <= z { [0.0, 1.0, 0.0] } else { [0.0, 0.0, 1.0] }
}
fn tool_flange_body(m, m6, pi) {
    let face = m6.g.len + 29.0;
    let body = foot_disc(m, m6);
    body = body.union(m.tool(|t| axial_disc(t, m6, m6.g.len + 22.0, face, 34.0, "spigot")),
        #{ label: "spigot" });
    body = body.subtract(m.tool(|t| axial_disc(t, m6, face - 6.0, face + 1.0, 15.75, "centring")),
        #{ label: "centring_bore" });
    body = iso_taps(body, m6, face, pi);
    let dowel = 3.0 * pi / 8.0;
    body = body.subtract(m.tool(|t| rod_at(t, m_point(m6.M, [25.0 * cos(dowel), 25.0 * sin(dowel), face + 1.0]),
        m_dir(m6.M, [0.0, 0.0, -1.0]), 3.0, 9.0, "dowel")), #{ label: "dowel" });
    body.subtract(m.tool(|t| cbore_tool(t, m6, m6.s.pcd / 2.0, m6.s.bolts, face, m6.cb_d + 6.0, 18.5, 10.0, pi)),
        #{ label: "flange_bolts" })
}

// ---- One joint's drive in its links ----
// The housing side (stator, windings, ring pins, outer race) in the
// link before the joint.
fn place_static(sub, p, md) {
    let z_axis = axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
    let j = "j" + (md.k + 1);
    place_at(sub.place(p.housing, j + "_housing"), md.M);
    place_at(sub.place(p.stator, j + "_stator"), md.M);
    for k in 0..12 {
        place_at(sub.place(p.winding, j + "_winding_" + k).rotate(z_axis, 360.0 * k / 12.0), md.M);
    }
    for k in 0..md.s.n {
        place_at(sub.place(p.ring_pin, j + "_ring_pin_" + k).rotate(z_axis, 360.0 * k / md.s.n), md.M);
    }
    place_at(sub.place(p.outer_race, j + "_outer_race"), md.M);
}
// The output side (flange, output pins, inner race) in the link after.
fn place_output(sub, p, md) {
    let z_axis = axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
    let j = "j" + (md.k + 1);
    place_at(sub.place(p.flange, j + "_flange"), md.M);
    for k in 0..6 {
        place_at(sub.place(p.output_pin, j + "_output_pin_" + k).rotate(z_axis, 60.0 * k), md.M);
    }
    place_at(sub.place(p.inner_race, j + "_inner_race"), md.M);
}
// The rotor and roller-cage sub-assemblies of one module size.
fn drive_defs(a, p, s, pi) {
    let g = p.g;
    let z_axis = axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
    let rotor = a.subassembly("rotor_" + s.name, |sub| {
        sub.place(p.shaft, "shaft");
        for k in 0..10 {
            sub.place(p.magnet, "magnet_" + k).rotate(z_axis, 36.0 * k);
        }
        sub.place(p.encoder, "encoder");
    });
    let rollers = ((2.0 * pi * g.rb) / (s.d + 1.2)).floor().to_int();
    let count = rollers - rollers % 2;
    let cage = a.subassembly("cage_" + s.name, |sub| {
        let tilt = axis_at([0.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        for k in 0..count {
            let sign = if k % 2 == 0 { 45.0 } else { -45.0 };
            sub.place(p.roller, "roller_" + k).rotate(tilt, sign).translate(g.rb, 0.0, g.zb)
                .rotate(z_axis, 360.0 * k / count);
        }
    });
    #{ rotor: rotor, cage: cage }
}
// Joint k's revolute between its static link st and output link out,
// and the drive parts that follow it.
fn place_joint(a, arm, md, p, defs, st, out) {
    let s = md.s;
    let g = md.g;
    let j = "j" + (md.k + 1);
    let z_axis = axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
    let ro = place_at(a.place(defs.rotor, j + "_rotor"), md.M);
    let da = place_at(a.place(p.disc, j + "_disc_a").translate(s.e, 0.0, g.z_d0), md.M);
    let db = place_at(a.place(p.disc, j + "_disc_b").rotate(z_axis, 180.0).translate(-s.e, 0.0, g.z_db), md.M);
    let cg = place_at(a.place(defs.cage, j + "_cage"), md.M);
    a.revolute(j, st.axis(md.o, md.z), out.axis(md.o, md.z),
        #{ zero: [st.datum("x"), out.datum("x")], angle: 0.0, limits: arm.limits[md.k] });
    // A revolute holds its two axis points level (axial offset 0): the
    // rotor and cage turn about the module's own origin.
    let mo = m_point(md.M, [0.0, 0.0, 0.0]);
    a.revolute(j + "_rotor", st.axis(mo, md.z), ro.axis([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
        #{ zero: [st.datum("x"), ro.datum("x")], follow: j, ratio: -(s.n - 1.0) });
    a.revolute(j + "_disc_a", ro.axis([s.e, 0.0, g.z_d0], [0.0, 0.0, 1.0]), da.datum("z"),
        #{ zero: [st.datum("x"), da.datum("x")], follow: j, ratio: 1.0 });
    a.revolute(j + "_disc_b", ro.axis([-s.e, 0.0, g.z_db], [0.0, 0.0, 1.0]), db.datum("z"),
        #{ zero: [st.datum("x"), db.datum("x")], follow: j, ratio: 1.0, phase: 180.0 });
    a.revolute(j + "_cage", st.axis(mo, md.z), cg.axis([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
        #{ zero: [st.datum("x"), cg.datum("x")], follow: j, ratio: 0.5 });
}
// A two-root link joined along its tube (x through the joint at
// (y = 0, z = zc)): the tube end slides to −param on the same line.
fn tube_mates(a, root, end, x_root, x_end, zc, param) {
    a.concentric(root.axis([x_root, 0.0, zc], [1.0, 0.0, 0.0]), end.axis([x_end, 0.0, zc], [1.0, 0.0, 0.0]));
    a.coincident(root.plane([x_root, 0.0, zc], [0.0, 0.0, 1.0]), end.plane([x_end, 0.0, zc], [0.0, 0.0, 1.0]),
        #{ aligned: true });
    a.distance(root.plane([x_root, 0.0, zc], [1.0, 0.0, 0.0]), end.plane([x_end, 0.0, zc], [1.0, 0.0, 0.0]),
        p(param), #{ aligned: true, scale: -1.0 });
}
// Fill a link: its castings; the output side of joint `prev` (−1:
// none) with the foot bolts; the housing side of joint `next` with its
// mount-ring bolts. The tool flange's bolts sit in its deeper pockets.
fn fill_link(sub, castings, prev, next, parts, mods, arm, bolts, pi) {
    for c in castings {
        sub.place(c[0], c[1]);
    }
    if prev >= 0 {
        let md = mods[prev];
        place_output(sub, parts[arm.sizes[prev]], md);
        let face = if prev == 5 { md.g.len + 29.0 } else { md.g.len + 23.0 };
        let depth = if prev == 5 { md.cb_d + 6.0 } else { md.cb_d };
        // Heads at face - depth; the flange's tapped face at len + 11.
        place_bolts(sub, bolts[md.screw], md, md.s.pcd / 2.0, md.s.bolts, face, depth,
            face - depth - (md.g.len + 11.0), "j" + (prev + 1) + "_foot", pi);
    }
    if next >= 0 {
        let md = mods[next];
        place_static(sub, parts[arm.sizes[next]], md);
        // Heads at len + 10 - cb_d; the housing lip's tapped face at len.
        place_bolts(sub, bolts[md.screw], md, md.rl, md.s.bolts, md.g.len + 10.0, md.cb_d,
            10.0 - md.cb_d, "j" + (next + 1) + "_mount", pi);
    }
}

// ---- Gripper (§4.2, PK-724): a parallel jaw on the tool flange ----
// World coordinates at home: the tool axis is the line x = tx, z = tz
// along −y; the flange face is y = ty. The body: an adapter disc on the
// flange (a spigot in its centring bore, four counterbored M6 screws
// into the ISO tapped holes, a dowel), then a jaw block with an open
// slide channel along x on its far face. Each finger: a carriage in the
// channel and a jaw beyond the block; finger B is finger A turned 180°
// about the tool axis.
fn gripper_layout(m6) {
    let face = m6.g.len + 29.0;
    let f = m_point(m6.M, [0.0, 0.0, face]);
    #{
        tx: f[0], ty: f[1], tz: f[2],
        face: face,
        adapter: 10.0, adapter_r: 34.0,
        block_y: 30.0, block_x: 62.0, block_z: 20.0,
        chan_z: 8.0, chan_d: 10.0,
        finger_w: 18.0, jaw_len: 60.0, jaw_z: 12.0, carriage_z: 7.5,
        grip0: 40.0, travel: 80.0,
    }
}
// The ISO 9409 screw positions on the module frame (PCD 50, 22.5° + 90k).
fn iso_angle(k, pi) { pi / 8.0 + pi * k / 2.0 }
fn gripper_body(m, m6, gl, pi) {
    let d = m_dir(m6.M, [0.0, 0.0, 1.0]);
    let a0 = gl.face;
    // Adapter disc, then the spigot into the centring bore (0.05 radial
    // clearance, 0.5 short of its floor).
    let body = axial_disc(m, m6, a0, a0 + gl.adapter, gl.adapter_r, "adapter");
    body = body.union(m.tool(|t| axial_disc(t, m6, a0 - 5.5, a0 + 1.0, 15.7, "spigot_sketch")),
        #{ label: "spigot" });
    // The jaw block, overlapping the adapter 0.5 (its faces never
    // coplanar with the disc's).
    let y_top = gl.ty - gl.adapter + 0.5;
    let y_bot = y_top - gl.block_y;
    body = body.union(m.tool(|t| box(t, gl.tx - gl.block_x, gl.tx + gl.block_x, y_bot, y_top,
        gl.tz - gl.block_z, gl.tz + gl.block_z, "block_sketch")), #{ label: "block" });
    // The slide channel, open along x and on the far face.
    body = body.subtract(m.tool(|t| box(t, gl.tx - gl.block_x - 1.0, gl.tx + gl.block_x + 1.0,
        y_bot - 1.0, y_bot + gl.chan_d, gl.tz - gl.chan_z, gl.tz + gl.chan_z, "channel_sketch")),
        #{ label: "channel" });
    // Four counterbored M6 holes on the ISO pattern, from the adapter's
    // far face toward the flange; the pockets open past the block side.
    let screws = m.tool(|t| {
        let all = ();
        for k in 0..4 {
            let ang = iso_angle(k, pi);
            let u = 25.0 * cos(ang);
            let v = 25.0 * sin(ang);
            let top = m_point(m6.M, [u, v, a0 + gl.adapter + 40.0]);
            let pocket = rod_at(t, top, m_dir(m6.M, [0.0, 0.0, -1.0]), 5.5, 46.5, "m6_pocket_" + k);
            let hole = rod_at(t, m_point(m6.M, [u, v, a0 + 4.0]), m_dir(m6.M, [0.0, 0.0, -1.0]), 3.3,
                5.0, "m6_clear_" + k);
            let both = pocket.union(hole);
            all = if k == 0 { both } else { all.union(both) };
        }
        all
    });
    body = body.subtract(screws, #{ label: "flange_screws" });
    // Dowel hole Ø6 matching the flange's.
    let dw = 3.0 * pi / 8.0;
    body.subtract(m.tool(|t| rod_at(t, m_point(m6.M, [25.0 * cos(dw), 25.0 * sin(dw), a0 - 1.0]),
        m_dir(m6.M, [0.0, 0.0, 1.0]), 3.0, 9.0, "dowel_hole")), #{ label: "dowel_hole" })
}
// Finger A at grip0: jaw face at x = tx + grip0/2, the finger reaching
// out to +x.
fn finger_body(m, gl) {
    let x0 = gl.tx + gl.grip0 / 2.0;
    let x1 = x0 + gl.finger_w;
    let y_block = gl.ty - gl.adapter + 0.5 - gl.block_y;
    let body = box(m, x0, x1, y_block - 0.5 - gl.jaw_len, y_block - 0.5, gl.tz - gl.jaw_z, gl.tz + gl.jaw_z,
        "jaw");
    body.union(m.tool(|t| box(t, x0 + 1.0, x1 - 1.0, y_block - 2.5, y_block + gl.chan_d - 0.5,
        gl.tz - gl.carriage_z, gl.tz + gl.carriage_z, "carriage_sketch")), #{ label: "carriage" })
}
// A finger slides along x on the body: its y and z planes held, its jaw
// face `sign` × grip/2 from the tool axis.
fn finger_mates(a, link6, finger, gl, sign) {
    // In the finger's own frame (finger B is turned 180° about the
    // tool axis, so its z and x normals face the other way in world).
    let at = [gl.tx, gl.ty - 20.0, gl.tz];
    a.coincident(link6.plane(at, [0.0, 1.0, 0.0]), finger.plane(at, [0.0, 1.0, 0.0]), #{ aligned: true });
    a.coincident(link6.plane(at, [0.0, 0.0, 1.0]), finger.plane(at, [0.0, 0.0, 1.0]), #{ aligned: sign > 0.0 });
    a.distance(link6.plane(at, [1.0, 0.0, 0.0]),
        finger.plane([gl.tx + gl.grip0 / 2.0, at[1], at[2]], [1.0, 0.0, 0.0]),
        p("grip"), #{ aligned: sign > 0.0, scale: 0.5 * sign });
}

// ---- Pose presets (§4.4) ----
fn home(a) {
    for k in 1..7 { a.set_param("j" + k, 0.0); }
}
fn ready(a) {
    let q = [0.0, -60.0, 100.0, -130.0, -90.0, 0.0];
    for k in 0..6 { a.set_param("j" + (k + 1), q[k]); }
}
fn reach(a) {
    let q = [45.0, -20.0, 30.0, -100.0, -90.0, 0.0];
    for k in 0..6 { a.set_param("j" + (k + 1), q[k]); }
}

// ---- The arm ----
let arm = arm_layout(pi);
let frames = dh_frames(arm);
let mods = joint_mods(arm, frames, sizes);
for s in sizes {
    if (s.n - 1) % 2 == 0 {
        throw "module " + s.name + ": " + s.n + " ring pins give an even lobe count; the "
            + "180° twin disc needs N − 1 odd";
    }
}
let parts = [];
let defs = [];
for s in sizes {
    let p = module_parts(a, s, pi);
    parts.push(p);
    defs.push(drive_defs(a, p, s, pi));
}
let zc = arm.d[0];
let ua = -arm.a[1];
let fa = -arm.a[2];
a.param("upper_arm", ua, "m");
a.param("forearm", fa, "m");

let base = a.part("base", |m| base_body(m, mods[0], pi));
let shoulder = a.part("shoulder", |m| shoulder_body(m, mods[0], mods[1], pi));
let upper_root = a.part("upper_arm_root", |m| arm_root_body(m, mods[1], arm.upper_y, arm.upper_r,
    60.0, -arm.upper_root, -48.0, pi));
// The tubes take their length from the arm parameters (the end
// castings stay fixed; the mates carry the far end).
let upper_tube = a.part("upper_arm_tube", |m| {
    m.param("upper_arm", ua, "m");
    m.param("tube_ends", arm.upper_root + arm.upper_end, "m");
    m.param("bore_ends", arm.upper_root + arm.upper_end - 2.0, "m");
    tube_body(m, -arm.upper_root, arm.upper_y, zc, arm.upper_r, arm.upper_ri,
        "upper_arm - tube_ends", "upper_arm - bore_ends");
});
let upper_end = a.part("upper_arm_end", |m| arm_end_body(m, mods[2], arm.upper_y, arm.upper_r,
    -ua + 60.0, -(ua - arm.upper_end), 30.0, pi));
let fore_root = a.part("forearm_root", |m| arm_root_body(m, mods[2], arm.fore_y, arm.fore_r,
    52.0, -ua - arm.fore_root, -ua - 40.5, pi));
let fore_tube = a.part("forearm_tube", |m| {
    m.param("forearm", fa, "m");
    m.param("tube_ends", arm.fore_root + arm.fore_end, "m");
    m.param("bore_ends", arm.fore_root + arm.fore_end - 2.0, "m");
    tube_body(m, -ua - arm.fore_root, arm.fore_y, zc, arm.fore_r, arm.fore_ri,
        "forearm - tube_ends", "forearm - bore_ends");
});
let fore_end = a.part("forearm_end", |m| arm_end_body(m, mods[3], arm.fore_y, arm.fore_r,
    -ua - fa + 51.0, -ua - (fa - arm.fore_end), 25.0, pi));
let wrist_1 = a.part("wrist_1", |m| wrist_1_body(m, mods[3], mods[4], pi));
let wrist_2 = a.part("wrist_2", |m| wrist_2_body(m, mods[4], mods[5], pi));
let tool_flange = a.part("tool_flange", |m| tool_flange_body(m, mods[5], pi));
let gl = gripper_layout(mods[5]);
let gripper = a.part("gripper_body", |m| gripper_body(m, mods[5], gl, pi));
let finger = a.part("finger", |m| finger_body(m, gl));
let dowel = a.part("dowel_pin", |m| rod_at(m, m_point(mods[5].M, [25.0 * cos(3.0 * pi / 8.0),
    25.0 * sin(3.0 * pi / 8.0), gl.face - 6.0]), m_dir(mods[5].M, [0.0, 0.0, 1.0]), 2.98, 12.0, "pin"));
a.param("grip", gl.grip0, "m");
a.limits("grip", 0.0, gl.travel);
a.drive("upper_arm", [upper_tube]);
a.drive("forearm", [fore_tube]);
// Each tube keeps at least 60 mm between its end castings.
a.limits("upper_arm", arm.upper_root + arm.upper_end + 60.0, 700.0);
a.limits("forearm", arm.fore_root + arm.fore_end + 60.0, 650.0);

let bolt = a.part("bolt", #{ size: "M5", length: 12.0, cosmetic: false }, |m, q| {
    m.bolt(q.size, q.length, #{ head: "socket", cosmetic_thread: q.cosmetic, label: "bolt" });
});
let bolts = #{ M5: bolt, M4: bolt.variant(#{ size: "M4", length: 12.0 }) };
let bolt_m8 = bolt.variant(#{ size: "M8", length: 20.0 });
let bolt_m6 = bolt.variant(#{ size: "M6", length: 12.0 });
for d in [base, shoulder, upper_root, upper_end, fore_root, fore_end, wrist_1, wrist_2] {
    a.appearance(d, #{ color: "white", metal: 0.1, rough: 0.5 });
}
for d in [upper_tube, fore_tube] {
    a.appearance(d, #{ color: "aluminium", metal: 0.7, rough: 0.3 });
}
a.appearance(tool_flange, #{ color: "steel", metal: 0.85, rough: 0.25 });
a.appearance(gripper, #{ color: "gray", metal: 0.3, rough: 0.5 });
a.appearance(finger, #{ color: "orange", metal: 0.1, rough: 0.6 });
a.appearance(dowel, #{ color: "steel", metal: 0.9, rough: 0.2 });
for d in [bolts.M5, bolts.M4, bolt_m8, bolt_m6] {
    a.appearance(d, #{ color: "black", metal: 0.6, rough: 0.4 });
}

// ---- Links: one rigid sub-assembly each (§4.2) ----
// Every link is authored at home in world coordinates and placed at
// identity: a link's placement IS its motion from home, and joint
// axes are the same literals in both links' frames.
let anchors = [[95.0, 95.0], [-95.0, 95.0], [-95.0, -95.0], [95.0, -95.0]];
let link_0 = a.subassembly("link_0", |sub| {
    fill_link(sub, [[base, "base"]], -1, 0, parts, mods, arm, bolts, pi);
    for k in 0..4 {
        let c = anchors[k];
        sub.place(bolt_m8, "anchor_" + k).frame([c[0], c[1], -74.5], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]);
    }
});
let link_1 = a.subassembly("link_1", |sub| {
    fill_link(sub, [[shoulder, "shoulder"]], 0, 1, parts, mods, arm, bolts, pi);
});
let link_2_tube = a.subassembly("link_2_tube", |sub| {
    fill_link(sub, [[upper_root, "root"], [upper_tube, "tube"]], 1, -1, parts, mods, arm, bolts, pi);
});
let link_2 = a.subassembly("link_2", |sub| {
    fill_link(sub, [[upper_end, "elbow"]], -1, 2, parts, mods, arm, bolts, pi);
});
let link_3_tube = a.subassembly("link_3_tube", |sub| {
    fill_link(sub, [[fore_root, "root"], [fore_tube, "tube"]], 2, -1, parts, mods, arm, bolts, pi);
});
let link_3 = a.subassembly("link_3", |sub| {
    fill_link(sub, [[fore_end, "wrist_end"]], -1, 3, parts, mods, arm, bolts, pi);
});
let link_4 = a.subassembly("link_4", |sub| {
    fill_link(sub, [[wrist_1, "wrist_1"]], 3, 4, parts, mods, arm, bolts, pi);
});
let link_5 = a.subassembly("link_5", |sub| {
    fill_link(sub, [[wrist_2, "wrist_2"]], 4, 5, parts, mods, arm, bolts, pi);
});
let link_6 = a.subassembly("link_6", |sub| {
    fill_link(sub, [[tool_flange, "tool_flange"], [gripper, "gripper"], [dowel, "dowel"]], 5, -1,
        parts, mods, arm, bolts, pi);
    // The gripper's four M6 screws, heads on their pocket floors 3.5
    // (3.5 leads) out of the ISO taps, which are cut in world
    // coordinates: half a turn from the taps' seed (PK-476, PK-876).
    let d = m_dir(mods[5].M, [0.0, 0.0, -1.0]);
    let s = insertion_seed(d);
    for k in 0..4 {
        let ang = iso_angle(k, pi);
        let at = m_point(mods[5].M, [25.0 * cos(ang), 25.0 * sin(ang), gl.face + 3.5]);
        sub.place(bolt_m6, "gripper_screw_" + k).frame(at, d, [-s[0], -s[1], -s[2]]);
    }
});

// ---- Joints ----
let roots = [];
for l in [[link_0, "link_0"], [link_1, "link_1"], [link_2_tube, "link_2_tube"], [link_2, "link_2"],
          [link_3_tube, "link_3_tube"], [link_3, "link_3"], [link_4, "link_4"], [link_5, "link_5"],
          [link_6, "link_6"]] {
    roots.push(a.place(l[0], l[1]));
}
a.ground(roots[0]);
// Joint k's static and output roots.
let st = [0, 1, 3, 5, 6, 7];
let out = [1, 2, 4, 6, 7, 8];
for k in 0..6 {
    let si = arm.sizes[k];
    place_joint(a, arm, mods[k], parts[si], defs[si], roots[st[k]], roots[out[k]]);
    if k == 1 { tube_mates(a, roots[2], roots[3], 0.0, -ua, zc, "upper_arm"); }
    if k == 2 { tube_mates(a, roots[4], roots[5], -ua, -ua - fa, zc, "forearm"); }
}
// Fingers: free roots sliding on the gripper body, opened by `grip`.
let finger_a = a.place(finger, "finger_a");
let finger_b = a.place(finger, "finger_b").rotate(axis_at([gl.tx, 0.0, gl.tz], [0.0, 1.0, 0.0]), 180.0);
finger_mates(a, roots[8], finger_a, gl, 1.0);
finger_mates(a, roots[8], finger_b, gl, -1.0);
a.assert_watertight();

// The displayed pose: home, ready or reach.
ready(a);
a
