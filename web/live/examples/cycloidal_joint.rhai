// Cycloidal joint modules for a six-axis arm, one of each size (L, M,
// S), as an ASSEMBLY (design docs/architecture/six-axis-arm.md §4.2,
// §5.2). Each module is a motor driving a cycloidal reducer: the rotor
// carries two eccentric cams 180° apart, each cam rolls a cycloid disc
// around N fixed ring pins, and six output pins through the discs turn
// the output flange at 1/(N−1) of the rotor speed. A cross-roller
// bearing carries the flange; its rollers sit at alternating ±45°.
//
// Every relation is driven by ONE slider per module (j_L, j_M, j_S): the
// output turns with it, the rotor follows at −(N−1), each disc follows
// 1:1 riding its cam's axis, and the roller cage follows at 1/2. Move a
// slider and the discs wobble round the pins; nothing rebuilds.
//
// The disc flank is the equidistant of an epitrochoid offset by the pin
// radius plus a backlash b. The sketch has no spline, so the curve is
// fitted with G1 biarcs, split until the sampled deviation is at most
// b/4: the flanks are exact cylinders and the backlash dominates the
// fit error. A backlash below the fit throws before anything is built.
//
// Fallback boundary (design §6, E1, PK-801): the bolt circles are
// tap-drill bores labelled `tapped_M5` / `tapped_M4`. A real helical
// thread costs ~1 s per hole (2026-10-14); a joint module has 24 and
// the arm ~150, past its 60 s cold-build budget. The
// encoder disc's 64 windows open at the rim (a slotted wheel): one
// sketch, no booleans.
let a = assembly("cycloidal_joint");
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
// `count` holes of radius hole_r on a circle of radius ring_r, from z
// for depth: one tool body.
fn bolt_tool(t, count, ring_r, hole_r, z, depth, pi) {
    let all = t.profile(xy_at(z)).circle(hole_r).at(ring_r, 0.0).sketch("tap_0").extrude(blind(depth));
    for k in 1..count {
        let ang = 2.0 * pi * k / count;
        let hole = t.profile(xy_at(z)).circle(hole_r).at(ring_r * cos(ang), ring_r * sin(ang))
            .sketch("tap_" + k).extrude(blind(depth));
        all = all.union(hole);
    }
    all
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
    // Bolt circle through the mounting lip: tap-drill bores (E1
    // fallback), one tool.
    let holes = m.tool(|t| bolt_tool(t, s.bolts, g.ro + 4.5, s.tap / 2.0, g.len - 9.0, 10.0, pi));
    body.subtract(holes, #{ label: "tapped_" + s.screw })
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
// with the output bolt circle (tap-drill bores, E1 fallback).
fn flange_body(m, g, s, pi) {
    let top = g.len + 1.0;
    let z0 = g.zb - g.half - 1.5;
    let p = m.profile(xz()).path(g.hub_ri, z0).line_to(g.hub_r, z0).line_to(g.hub_r, top);
    p = p.line_to(g.flange_r, top).line_to(g.flange_r, top + 10.0).line_to(0.0, top + 10.0);
    p = p.line_to(0.0, top).line_to(g.hub_ri, top).close();
    let body = p.sketch("flange").revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle());
    let holes = m.tool(|t| bolt_tool(t, s.bolts, s.pcd / 2.0, s.tap / 2.0, top + 2.0, 9.0, pi));
    body.subtract(holes, #{ label: "tapped_" + s.screw })
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

// ---- One module: sub-assemblies, placements and joints ----
// The module sits at x = at; its slider is "j_" + size.
fn place_module(a, s, p, at, pi) {
    let g = p.g;
    let z_axis = axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
    let static_side = a.subassembly("static_" + s.name, |sub| {
        sub.place(p.housing, "housing");
        sub.place(p.stator, "stator");
        for k in 0..12 {
            sub.place(p.winding, "winding_" + k).rotate(z_axis, 360.0 * k / 12.0);
        }
        for k in 0..s.n {
            sub.place(p.ring_pin, "ring_pin_" + k).rotate(z_axis, 360.0 * k / s.n);
        }
        sub.place(p.outer_race, "outer_race");
    });
    let output_side = a.subassembly("output_side_" + s.name, |sub| {
        sub.place(p.flange, "flange");
        for k in 0..6 {
            sub.place(p.output_pin, "pin_" + k).rotate(z_axis, 60.0 * k);
        }
        sub.place(p.inner_race, "inner_race");
    });
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

    let st = a.place(static_side, "static_" + s.name).translate(at, 0.0, 0.0);
    let out = a.place(output_side, "output_" + s.name).translate(at, 0.0, 0.0);
    let ro = a.place(rotor, "rotor_" + s.name).translate(at, 0.0, 0.0);
    let da = a.place(p.disc, "disc_" + s.name + "_a").translate(at + s.e, 0.0, g.z_d0);
    // Disc B is the same disc turned 180° on the opposite cam: with an
    // odd lobe count (N − 1) that is half a lobe for the flank, and the
    // six windows map onto themselves.
    let db = a.place(p.disc, "disc_" + s.name + "_b")
        .rotate(z_axis, 180.0).translate(at - s.e, 0.0, g.z_db);
    let cg = a.place(cage, "cage_" + s.name).translate(at, 0.0, 0.0);
    a.ground(st);

    // One slider per module; everything else follows it.
    let j = "j_" + s.name;
    a.revolute(j, st.axis([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), out.axis([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
        #{ zero: [st.datum("x"), out.datum("x")], angle: 0.0, limits: [-360.0, 360.0] });
    a.revolute(j + "_rotor", st.axis([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), ro.axis([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
        #{ zero: [st.datum("x"), ro.datum("x")], follow: j, ratio: -(s.n - 1.0) });
    a.revolute(j + "_disc_a", ro.axis([s.e, 0.0, g.z_d0], [0.0, 0.0, 1.0]), da.datum("z"),
        #{ zero: [st.datum("x"), da.datum("x")], follow: j, ratio: 1.0 });
    a.revolute(j + "_disc_b", ro.axis([-s.e, 0.0, g.z_db], [0.0, 0.0, 1.0]), db.datum("z"),
        #{ zero: [st.datum("x"), db.datum("x")], follow: j, ratio: 1.0, phase: 180.0 });
    a.revolute(j + "_cage", st.axis([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), cg.axis([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
        #{ zero: [st.datum("x"), cg.datum("x")], follow: j, ratio: 0.5 });
}

let x = 0.0;
for s in sizes {
    if (s.n - 1) % 2 == 0 {
        throw "module " + s.name + ": " + s.n + " ring pins give an even lobe count; the "
            + "180° twin disc needs N − 1 odd";
    }
    let parts = module_parts(a, s, pi);
    place_module(a, s, parts, x, pi);
    x += s.od + 40.0;
}
a
