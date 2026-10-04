// A centrifugal impeller, built the way one is machined: a solid
// blank, then seven spiral flow channels carved out between the
// blades. The channel walls are logarithmic spirals sampled by this
// script — the blade shape IS the leftover metal.
let m = model();
m.param("bore_r", 0.1, "m");

let pi = 3.141592653589793;

// ---- Blank: disc, centre boss, bore ----
let wheel = m.profile(xy()).circle(1.0).sketch("blank_disc").extrude(blind(0.42));
let boss = m.tool(|m| m.profile(xy_at(0.4)).circle(0.22).sketch("hub_boss").extrude(blind(0.15)));
wheel = wheel.union(boss);
let bore = m.tool(|m| m.profile(xy_at(-0.1)).circle(p("bore_r")).sketch("bore_sketch").extrude(blind(0.75)));
wheel = wheel.subtract(bore);

// ---- The flow cavity: one connected star-shaped cut ----
// Seven channel arms joined by root caps around each blade — a
// single computed polygon, one subtract. Each channel wall follows
// the logarithmic spiral r(t) = r0 * (r1/r0)^t wound back by `wrap`
// radians; the blade between two walls tapers from a thick root to
// a near-sharp trailing edge, the way a pump wheel's blades do.
let blades = 7;
let r0 = 0.30;
let r1 = 0.95;
let wrap = 1.1;
let half_t = 0.05;
let steps = 16;
let pitch = 2.0 * pi / blades;
// Walls run to r = 0.97, just inside the rim: the wheel keeps its
// shroud ring.
let t_max = ln(0.97 / r0) / ln(r1 / r0);

// Blade half-thickness, tapered toward the trailing edge.
fn ht(t, t_max, half_t) {
    let f = (t_max - t) / 0.35;
    let f = if f > 1.0 { 1.0 } else if f < 0.0 { 0.0 } else { f };
    half_t * (0.24 + 0.76 * f)
}

fn wall_point(base, t, side, r0, r1, wrap, half_t, t_max) {
    let r = r0 * (r1 / r0) ** t;
    let theta = base + wrap * t;
    let dr = r * ln(r1 / r0);
    let dtheta = wrap * r;
    let len = sqrt(dr * dr + dtheta * dtheta);
    let h = ht(t, t_max, half_t);
    let nr = -dtheta / len * side * h;
    let nt = dr / len * side * h;
    [(r + nr) * cos(theta) - nt * sin(theta),
     (r + nr) * sin(theta) + nt * cos(theta)]
}

// A point f of the way from `from` to `to`, interpolating radius and
// angle separately (the short way round) — closure arcs that stay
// clear of the walls they join.
fn polar_blend(from, to, f, pi) {
    let r_a = sqrt(from[0] * from[0] + from[1] * from[1]);
    let r_b = sqrt(to[0] * to[0] + to[1] * to[1]);
    let a_a = atan(from[1] / from[0]) + if from[0] < 0.0 { pi } else { 0.0 };
    let a_b = atan(to[1] / to[0]) + if to[0] < 0.0 { pi } else { 0.0 };
    let a_b = if a_b - a_a > pi { a_b - 2.0 * pi }
              else if a_a - a_b > pi { a_b + 2.0 * pi }
              else { a_b };
    let r = r_a + (r_b - r_a) * f;
    let a = a_a + (a_b - a_a) * f;
    [r * cos(a), r * sin(a)]
}

let outline = [];
for b in 0..blades {
    // Out along one wall of the channel...
    for s in 0..=steps {
        let t = t_max * s / steps;
        outline.push(wall_point(b * pitch, t, 1.0, r0, r1, wrap, half_t, t_max));
    }
    // ...around the channel mouth...
    let from = wall_point(b * pitch, t_max, 1.0, r0, r1, wrap, half_t, t_max);
    let to = wall_point(b * pitch + pitch, t_max, -1.0, r0, r1, wrap, half_t, t_max);
    for s in 1..6 {
        outline.push(polar_blend(from, to, s * 1.0 / 6.0, pi));
    }
    // ...back along the next blade's wall...
    for s in 0..=steps {
        let t = t_max * (steps - s) / steps;
        outline.push(wall_point(b * pitch + pitch, t, -1.0, r0, r1, wrap, half_t, t_max));
    }
    // ...and around that blade's root to the next arm.
    let from = wall_point(b * pitch + pitch, 0.0, -1.0, r0, r1, wrap, half_t, t_max);
    let to = wall_point(b * pitch + pitch, 0.0, 1.0, r0, r1, wrap, half_t, t_max);
    for s in 1..4 {
        outline.push(polar_blend(from, to, s * 1.0 / 4.0, pi));
    }
}
let cavity = m.tool(|m| m.profile(xy_at(0.12)).polygon(outline)
    .sketch("flow_cavity")
    .extrude(blind(0.5)));
wheel = wheel.subtract(cavity);

m.assert_watertight();
m
