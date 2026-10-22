//! PK-459 stage 2: booleans on a doubly-wrapped torus face serve
//! watertight with their true volume
//! (docs/architecture/torus-cut-domain-chart-walk.md §10).
//!
//! A circle revolved a full circle is ONE torus face bounded by its
//! seam square (a u-seam meridian and a v-seam latitude through one
//! vertex). The chart walk takes it on its cut domain, the seams as
//! the domain's borders. Measured before: every radial drill declined
//! or failed validation, and a drill clear of both seams served
//! watertight with volume 0.000209 (the torus mis-classified away).
//!
//! The cases cover each way a record meets the domain: through the
//! seam corner, across the v-seam, across the u-seam, clear of both
//! (island holes), non-separating meridian wraps (a half torus),
//! latitude wraps (a slab cut) and a union. Oracles are exact
//! integrals; the 1% gate absorbs the mesher's chord bias (0.16% at
//! this chord height).

use popkern_rhai::ScriptHost;

const R: f64 = 0.4;
const RT: f64 = 0.1;
const RD: f64 = 0.02;

const TORUS: &str = r#"let m = model();
let b = m.profile(xz()).circle(0.1).at(0.4, 0.0).sketch("ring")
    .revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle());
"#;

/// The torus with a radial drill of radius `RD` entering from azimuth
/// `theta` at height `z0`, through both tube walls.
fn drill_script(theta: f64, z0: f64) -> String {
    let (c, s) = (theta.cos(), theta.sin());
    format!(
        "{TORUS}b = b.subtract(m.tool(|t| t.profile(plane_at([{x}, {y}, {z0}], [{nx}, {ny}, 0.0])).circle({RD}).sketch(\"bore\").extrude(blind(0.4))), #{{ label: \"drill\" }});\nm\n",
        x = 0.6 * c,
        y = 0.6 * s,
        nx = -c,
        ny = -s
    )
}

/// Mesh volume and open-segment census at the 1e-9 quantum.
fn census(script: &str) -> (f64, usize) {
    let mut host = ScriptHost::new();
    let outcome = host
        .run(script, "torus.rhai")
        .unwrap_or_else(|e| panic!("{e}\n{}", e.payload));
    let mut model = outcome.model.expect("model");
    let part = model.to_part().expect("part");
    let brep = part.brep();
    brep.validate().expect("valid result");
    let options = popkern_mesh::TessellationOptions {
        chord_height: 5.0e-4,
        min_subdivisions: 4,
        max_subdivisions: 512,
        ..Default::default()
    };
    let mesh = popkern_mesh::tessellate(brep, &options).expect("mesh");
    let mut vol6 = 0.0;
    type QuantizedPoint = (i64, i64, i64);
    let mut uses: std::collections::HashMap<(QuantizedPoint, QuantizedPoint), i64> =
        Default::default();
    let key = |p: popkern_core::geometry::Point3| {
        (
            (p.x * 1e9).round() as i64,
            (p.y * 1e9).round() as i64,
            (p.z * 1e9).round() as i64,
        )
    };
    for t in &mesh.triangles {
        let [a, b, c] = t.map(|i| mesh.vertices[i].position);
        vol6 += a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
            + a.z * (b.x * c.y - b.y * c.x);
        for (i, j) in [(0usize, 1usize), (1, 2), (2, 0)] {
            let a = key(mesh.vertices[t[i]].position);
            let b = key(mesh.vertices[t[j]].position);
            *uses.entry((a.min(b), a.max(b))).or_insert(0) += if a < b { 1 } else { -1 };
        }
    }
    let open = uses.values().filter(|&&c| c != 0).count();
    (vol6.abs() / 6.0, open)
}

fn torus_volume() -> f64 {
    2.0 * std::f64::consts::PI.powi(2) * R * RT * RT
}

/// The torus minus a radial drill at height `z0` (exact: the drill's
/// chord through the tube, integrated over its disc).
fn drilled_volume(z0: f64) -> f64 {
    let n = 400;
    let h = 2.0 * RD / n as f64;
    let mut removed = 0.0;
    for i in 0..n {
        for j in 0..n {
            let (s, z) = (-RD + (i as f64 + 0.5) * h, -RD + (j as f64 + 0.5) * h);
            if s * s + z * z >= RD * RD {
                continue;
            }
            let w = (RT * RT - (z + z0).powi(2)).sqrt();
            removed +=
                (((R + w).powi(2) - s * s).sqrt() - ((R - w).powi(2) - s * s).sqrt()) * h * h;
        }
    }
    torus_volume() - removed
}

/// Area of the disc of radius `a` beyond the chord at distance `c`.
fn disc_segment(a: f64, c: f64) -> f64 {
    if a <= c {
        0.0
    } else {
        a * a * (c / a).acos() - c * (a * a - c * c).sqrt()
    }
}

fn assert_serves(script: &str, want: f64, label: &str) {
    let (volume, open) = census(script);
    assert_eq!(open, 0, "[{label}] torn mesh: {open} open segments");
    let err = (volume - want).abs() / want;
    assert!(
        err < 0.01,
        "[{label}] volume {volume:.6} vs exact {want:.6} ({:.2}%)",
        err * 100.0
    );
}

#[test]
fn radial_drills_serve_at_seam_and_generic_stations() {
    // (azimuth, height): the seam corner, the v-seam (the outer
    // equator) at two azimuths, the u-seam alone, and clear of both.
    for (theta, z0) in [(0.0, 0.0), (0.7, 0.0), (2.1, 0.0), (0.0, 0.05), (0.7, 0.05)] {
        assert_serves(
            &drill_script(theta, z0),
            drilled_volume(z0),
            &format!("drill u {theta} z {z0}"),
        );
    }
}

#[test]
fn meridian_wraps_serve_a_half_torus() {
    // x > 0: two meridian records, each a non-separating wrap; the kept
    // half straddles the u-seam.
    let script = format!(
        "{TORUS}b = b.intersect(m.tool(|t| t.profile(xy_at(-0.5)).rect(1.0, 2.0).at(0.5, 0.0).sketch(\"bx\").extrude(blind(1.0))));\nm\n"
    );
    assert_serves(&script, 0.5 * torus_volume(), "half torus");
}

#[test]
fn latitude_and_chord_cuts_serve() {
    // z > 0.03 removed: two latitude records wrapping u.
    let d = 0.03;
    let slab = format!(
        "{TORUS}b = b.subtract(m.tool(|t| t.profile(xy_at({d})).rect(2.0, 2.0).sketch(\"bx\").extrude(blind(1.0))));\nm\n"
    );
    let kept = std::f64::consts::PI * RT * RT - disc_segment(RT, d);
    assert_serves(&slab, 2.0 * std::f64::consts::PI * R * kept, "slab");

    // x > 0.4 removed: the cut plane meets the tube around the u-seam.
    let cut = format!(
        "{TORUS}b = b.subtract(m.tool(|t| t.profile(xy_at(-0.5)).rect(1.0, 2.0).at(0.9, 0.0).sketch(\"bx\").extrude(blind(1.0))));\nm\n"
    );
    let n = 20000;
    let h = 2.0 * RT / n as f64;
    let mut want = 0.0;
    for i in 0..n {
        let z = -RT + (i as f64 + 0.5) * h;
        let w = (RT * RT - z * z).sqrt();
        let annulus = std::f64::consts::PI * ((R + w).powi(2) - (R - w).powi(2));
        want += (annulus - (disc_segment(R + w, 0.4) - disc_segment(R - w, 0.4))) * h;
    }
    assert_serves(&cut, want, "chord cut");
}

#[test]
fn a_cube_unions_across_the_v_seam() {
    // A 0.1 cube centred on the outer equator (the v-seam) at u = 0.
    let script = format!(
        "{TORUS}b = b.union(m.tool(|t| t.profile(xy_at(-0.05)).rect(0.1, 0.1).at(0.5, 0.0).sketch(\"bx\").extrude(blind(0.1))));\nm\n"
    );
    let n = 1000;
    let h = 0.1 / n as f64;
    let mut overlap = 0.0;
    for i in 0..n {
        let y = -0.05 + (i as f64 + 0.5) * h;
        for j in 0..n {
            let z = -0.05 + (j as f64 + 0.5) * h;
            let w = (RT * RT - z * z).sqrt();
            let lo = ((R - w).powi(2) - y * y).max(0.0).sqrt().max(0.45);
            let hi = ((R + w).powi(2) - y * y).sqrt().min(0.55);
            if hi > lo {
                overlap += (hi - lo) * h * h;
            }
        }
    }
    assert_serves(&script, torus_volume() + 0.001 - overlap, "cube union");
}
