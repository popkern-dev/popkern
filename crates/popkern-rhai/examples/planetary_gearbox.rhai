// A 4:1 planetary gearbox: a keyed sun gear, three planets, an internal
// ring gear, and a three-arm spider carrier with an output boss — six
// meshing bodies from one script. Every tooth flank is a TRUE involute
// computed below, and one outline function cuts both the external gears
// and the internal ring. The planets are PHASED so all three mesh
// simultaneously; the assembly constraint (Zs + Zr) % N == 0 is checked
// in the source.
let m = model();
m.param("width", 0.4, "m");
m.param("pin_bore", 0.075, "m");

let pi = 3.141592653589793;
let pressure = 30.0 * pi / 180.0; // 30°: undercut-free down to 8 teeth
let z_sun = 9;
let z_planet = 9;
let planets = 3;
let z_ring = z_sun + 2 * z_planet; // 27 — concentric by construction
// Ratio with the ring held, sun driving, carrier out: 1 + Zr/Zs = 4.0
if (z_sun + z_ring) % planets != 0 {
    throw "planet assembly constraint violated: (Zs + Zr) % N != 0";
}
let tooth_module = 1.5 / z_ring; // ring pitch radius 0.75
let backlash = 0.004; // angular flank thinning per flank

fn inv(a) { tan(a) - a }
// Involute roll angle at radius r off base circle rb. At r == rb the
// ratio rounds a hair past 1 (and the roll angle is zero there anyway),
// so return zero rather than letting acos go NaN — note Rhai's `>` is
// approximate, so guard with `>=` which counts the rounded-equal case.
fn roll(rb, r) { let x = rb / r; if x >= 1.0 { 0.0 } else { acos(x) } }

// Involute-toothed outline, tooth centers at phase + k*2pi/z: flanks
// run from max(base, root) radius to the tip in `flank_steps` chords
// (the external gears are cheap and take more; the ring cutter drives
// the boolean cost, so it takes fewer) — starting below the root
// circle would notch a slit into every valley — and the root is
// walked at rf. widen < 0 thins teeth (backlash on the external
// gears); widen > 0 widens them (the ring cutter's teeth are the
// ring's tooth SPACES, so widening the cutter is ring backlash).
// (Arithmetic is broken into small locals throughout: the scripting
// engine caps per-expression AST depth, and deep `a + b*c/d` chains
// nested in loops exceed it.)
fn toothed_outline(z, rb, ra, rf, cx, cy, widen, phase, pressure, pi, flank_steps) {
    let inv_p = inv(pressure);
    let half_p = pi / (2.0 * z);
    let pitch = 2.0 * pi / z;
    let r_lo = if rb > rf { rb } else { rf };
    let a_lo = roll(rb, r_lo);
    let half_lo = half_p + inv_p - inv(a_lo) + widen;
    // The valley is a real gap only when the flank starts above the
    // root circle; for the ring cutter (r_lo == rf) it collapses to a
    // single chord, so emit no interior valley points there.
    let has_valley = r_lo - rf >= 1e-12;
    let v_last = if has_valley { 3 } else { -1 };
    let pts = [];
    for k in 0..z {
        let rot = phase + k * pitch;
        // Right flank, root to tip.
        for s in 0..=flank_steps {
            let r = r_lo + (ra - r_lo) * s / flank_steps;
            let a = roll(rb, r);
            let half = half_p + inv_p - inv(a) + widen;
            let ang = rot - half;
            let qx = cx + r * cos(ang);
            let qy = cy + r * sin(ang);
            pts.push([qx, qy]);
        }
        // Left flank, tip back to root.
        for s in 0..=flank_steps {
            let r = ra - (ra - r_lo) * s / flank_steps;
            let a = roll(rb, r);
            let half = half_p + inv_p - inv(a) + widen;
            let ang = rot + half;
            let qx = cx + r * cos(ang);
            let qy = cy + r * sin(ang);
            pts.push([qx, qy]);
        }
        // Root valley walked to the next tooth.
        let vf = rot + half_lo;
        let vt = rot + pitch - half_lo;
        let vspan = vt - vf;
        for s in 0..=v_last {
            let th = vf + vspan * s / 3;
            let qx = cx + rf * cos(th);
            let qy = cy + rf * sin(th);
            pts.push([qx, qy]);
        }
    }
    pts
}

// ---- Sun gear (input, keyed bore) ----
let rp_s = tooth_module * z_sun / 2.0;
let rb_s = rp_s * cos(pressure);
let ra_s = rp_s + tooth_module;
let rf_s = rp_s - 1.25 * tooth_module;
let sun_outline = toothed_outline(z_sun, rb_s, ra_s, rf_s, 0.0, 0.0,
    -backlash, 0.0, pressure, pi, 6);
let sun = m.profile(xy()).polygon(sun_outline).sketch("sun_gear")
    .extrude(blind(p("width")));
let sun_bore = m.tool(|m| {
    // Bore and keyway are ONE tool with staggered z-extents.
    let hole = m.profile(xy_at(-0.1)).circle(0.1).sketch("sun_bore").extrude(blind(0.6));
    let key = m.profile(xy_at(-0.14)).rect(0.06, 0.05).at(0.0, 0.115)
        .sketch("sun_key").extrude(blind(0.68));
    hole.union(key)
});
sun = sun.subtract(sun_bore);

// ---- Three planets, phased to mesh ----
// A planet at carrier angle phi meshes the sun when its rotation is
// psi = phi*(1 + Zs/Zp) + pi + pi/Zp — tooth centers of one land on
// gap centers of the other along the line of centers.
let carrier_r = rp_s + tooth_module * z_planet / 2.0;
for k in 0..planets {
    let phi = k * 2.0 * pi / planets;
    let psi = phi * (1.0 + z_sun * 1.0 / z_planet) + pi + pi / z_planet;
    let cx = carrier_r * cos(phi);
    let cy = carrier_r * sin(phi);
    let outline = toothed_outline(z_planet, rb_s, ra_s, rf_s, cx, cy,
        -backlash, psi, pressure, pi, 6);
    let planet = m.profile(xy()).polygon(outline).sketch("planet_" + k)
        .extrude(blind(p("width")));
    let pin_hole = m.tool(|m| m.profile(xy_at(-0.1)).circle(p("pin_bore")).at(cx, cy)
        .sketch("planet_bore_" + k).extrude(blind(0.6)));
    planet = planet.subtract(pin_hole);
}

// ---- Ring gear: internal teeth cut by a toothed disc ----
// The cutter's teeth are the ring's tooth spaces: involutes of the
// ring base circle, tips at the ring root, valleys at the ring tooth
// tips. Its phase follows planet 0's: the (Zs+Zr)%N check above is
// exactly what lets the other two planets mesh the same ring.
let rp_r = tooth_module * z_ring / 2.0;
let rb_r = rp_r * cos(pressure);
let ring = m.profile(xy()).circle(0.92).sketch("ring_blank")
    .extrude(blind(p("width")));
let ring_phase = (pi + pi / z_planet) * z_planet / z_ring;
// The ring cutter drives the boolean cost, so it takes fewer flank
// chords (its teeth are small and peripheral); coarser chords fall
// INSIDE the true involute, thinning the ring teeth and only widening
// the mesh clearance.
let cutter_outline = toothed_outline(z_ring, rb_r, rp_r + 1.25 * tooth_module,
    rp_r - tooth_module, 0.0, 0.0, backlash, ring_phase, pressure, pi, 3);
let cutter = m.tool(|m| m.profile(xy_at(-0.1)).polygon(cutter_outline)
    .sketch("ring_cutter").extrude(blind(0.6)));
ring = ring.subtract(cutter);

// ---- Carrier: a three-arm spider so the gear train shows through ----
// A hub over the sun with three radial arms out to the planet pins,
// built as ONE polygon (arm-into-hub unions hit the split-precondition
// class when a straight arm crosses the round hub). The open spaces
// between the arms are what let the meshing gears be seen.
let hub_z = 0.44;
let hub_r = 0.24;
let arm_hw = 0.09; // arm half-width
let arm_out = carrier_r + 0.02; // reach just past the pins
// Each arm meets the hub circle where the offset edge is at radius
// hub_r: t0 = sqrt(hub_r^2 - hw^2), at hub angle +/- da.
let t0 = sqrt(hub_r * hub_r - arm_hw * arm_hw);
let da = atan(arm_hw / t0);
let spider = [];
for k in 0..planets {
    let phi = k * 2.0 * pi / planets;
    let ux = cos(phi);
    let uy = sin(phi);
    let px = -uy; // unit perpendicular
    let py = ux;
    let wx = arm_hw * px;
    let wy = arm_hw * py;
    // Out along the right edge, across the tip, back along the left.
    spider.push([t0 * ux - wx, t0 * uy - wy]);
    spider.push([arm_out * ux - wx, arm_out * uy - wy]);
    spider.push([arm_out * ux + wx, arm_out * uy + wy]);
    spider.push([t0 * ux + wx, t0 * uy + wy]);
    // Hub arc to the next arm.
    let a0 = phi + da;
    let a1 = (k + 1) * 2.0 * pi / planets - da;
    for s in 1..6 {
        let a = a0 + (a1 - a0) * s / 6;
        spider.push([hub_r * cos(a), hub_r * sin(a)]);
    }
}
let carrier = m.profile(xy_at(hub_z)).polygon(spider).sketch("carrier_spider")
    .extrude(blind(0.12));
for k in 0..planets {
    let phi = k * 2.0 * pi / planets;
    let pin = m.tool(|m| m.profile(xy_at(0.02)).circle(0.07)
        .at(carrier_r * cos(phi), carrier_r * sin(phi))
        .sketch("pin_" + k).extrude(blind(0.48)));
    carrier = carrier.union(pin);
}
let boss = m.tool(|m| m.profile(xy_at(0.52)).circle(0.16)
    .sketch("output_boss").extrude(blind(0.28)));
carrier = carrier.union(boss);
// A plain coaxial output bore through the boss and hub. (The keyed
// coupling lives on the SUN input; a keyway cut through this
// hub/boss union junction reaches the open chain-into-junction split
// class — PK-428 — so the carrier takes a plain bore.)
let out_bore = m.tool(|m| m.profile(xy_at(0.40)).circle(0.1)
    .sketch("output_bore").extrude(blind(0.45)));
carrier = carrier.subtract(out_bore);

m.assert_watertight();
m
