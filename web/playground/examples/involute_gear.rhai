// A spur gear whose tooth flanks are TRUE involutes, computed
// point-by-point in this script — the gear math lives in the source,
// not in a library. Twenty teeth, tooth_module 0.05, 20° pressure angle;
// bored, keyed, and chamfered like a real part.
let m = model();
m.param("width", 0.3, "m");

let pi = 3.141592653589793;
let teeth = 20;
let tooth_module = 0.05;
let pressure = 20.0 * pi / 180.0;

let rp = tooth_module * teeth / 2.0;        // pitch radius
let rb = rp * cos(pressure);          // base circle
let ra = rp + tooth_module;                 // addendum (tip)
let rf = rp - 1.25 * tooth_module;          // dedendum (root)

// Involute rollout angle and polar half-thickness of a flank point
// at radius r: at the pitch circle the half-thickness is pi/(2z),
// and it shrinks toward the tip by the involute function
// inv(a) = tan(a) - a.
fn inv(a) { tan(a) - a }
let inv_p = inv(pressure);
let half_p = pi / (2.0 * teeth);

// One tooth's outline points (built at angle 0, CCW: right flank
// root->tip, tip arc, left flank tip->root), then rotated per tooth.
let flank_steps = 8;
let outline = [];
for k in 0..teeth {
    let rot = k * 2.0 * pi / teeth;
    // Right flank: root to tip.
    for s in 0..=flank_steps {
        let r = rb + (ra - rb) * s / flank_steps;
        let alpha = acos(rb / r);
        let half = half_p + inv_p - inv(alpha);
        let theta = rot - half;
        outline.push([r * cos(theta), r * sin(theta)]);
    }
    // Left flank: tip back to root.
    for s in 0..=flank_steps {
        let r = ra - (ra - rb) * s / flank_steps;
        let alpha = acos(rb / r);
        let half = half_p + inv_p - inv(alpha);
        let theta = rot + half;
        outline.push([r * cos(theta), r * sin(theta)]);
    }
    // Root valley: drop to the dedendum circle and walk it to the
    // next tooth's right flank.
    let valley_from = rot + half_p + inv_p;
    let valley_to = rot + 2.0 * pi / teeth - (half_p + inv_p);
    for s in 0..=3 {
        let theta = valley_from + (valley_to - valley_from) * s / 3;
        outline.push([rf * cos(theta), rf * sin(theta)]);
    }
}

let gear = m.profile(xy()).polygon(outline).sketch("gear_outline")
    .extrude(blind(p("width")));

// The bore and its keyway are ONE tool — their union is a single
// closed outline, cut in one subtract.
let bore = m.tool(|m| {
    let hole = m.profile(xy_at(-0.1)).circle(0.15).sketch("bore_sketch").extrude(blind(0.5));
    // Staggered z-extents: both pierce the gear, and offset caps keep
    // the tool union transversal.
    let key = m.profile(xy_at(-0.14)).rect(0.08, 0.06).at(0.0, 0.17).sketch("key_sketch").extrude(blind(0.58));
    hole.union(key)
});
gear = gear.subtract(bore);

m.assert_watertight();
m
