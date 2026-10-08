//! Perpendicular cylinder tees (PK-727, PK-728): a main cylinder
//! r = 40 along z unioned with a branch along x at z = 60, at
//! millimetre scale — the castings of the six-axis arm.
//!
//! - PK-728: the analytic cylinder×cylinder curve was fitted from a
//!   fixed 512 samples; at this scale it strayed 2.4e-6 from the exact
//!   seam crossing (past the 2e-6 node test), the record was never
//!   split there and the main face tore (1118 open mesh edges). The
//!   fit is now densified to a tenth of the linear tolerance.
//! - PK-727: equal radii make the intersection two ellipses CROSSING at
//!   two singular points; the scan joined them into one self-touching
//!   loop the clip dropped, and the union SERVED both operands
//!   unmerged (two shells, volume = their sum). The branches now split
//!   at the crossings; a configuration the walk cannot yet build
//!   declines (PK-732) — never a wrong serve.
//!
//! Volumes are checked against the analytic value: the overlap of a
//! branch of radius `rb` through the main cylinder is
//! `∫ 2√(R²−y²) · 2√(rb²−y²) dy` over `|y| < rb`.

use popkern_rhai::ScriptHost;

const R: f64 = 40.0;
const H: f64 = 120.0;

fn tee(rb: f64, x0: f64, len: f64) -> String {
    format!(
        r#"
let m = model();
let main = m.profile(xy()).circle({R:?}).sketch("main").extrude(blind({H:?}));
let branch = m.tool(|t| t.profile(yz_at({x0:?})).circle({rb:?}).at(0.0, 60.0).sketch("branch").extrude(blind({len:?})));
let b = main.union(branch, #{{ label: "tee" }});
m.assert_watertight();
m
"#
    )
}

/// The volume of main ∪ branch for a branch spanning `x ∈ [x0, x0 + len]`
/// with `x0 ≤ −R` (a full cross) or `x0` inside the main cylinder such
/// that the branch end disk lies inside it.
fn analytic(rb: f64, x0: f64, len: f64) -> f64 {
    let n = 200_000;
    let mut overlap = 0.0;
    for i in 0..n {
        let y = -rb + 2.0 * rb * (f64::from(i) + 0.5) / f64::from(n);
        let half_chord = (R * R - y * y).sqrt();
        let depth = if x0 <= -R {
            2.0 * half_chord
        } else {
            half_chord - x0
        };
        overlap += depth * 2.0 * (rb * rb - y * y).sqrt() * 2.0 * rb / f64::from(n);
    }
    std::f64::consts::PI * (R * R * H + rb * rb * len) - overlap
}

/// Run the tee: `Some(volume)` when it serves (one watertight body),
/// `None` when the kernel declines.
fn served_volume(rb: f64, x0: f64, len: f64) -> Option<f64> {
    let mut host = ScriptHost::new();
    let outcome = host.run(&tee(rb, x0, len), "tee.rhai").ok()?;
    let mut model = outcome.model.expect("model");
    let part = model.to_part().expect("part");
    let brep = part.brep();
    Some(
        popkern_measure::body_volume(brep, &brep.tolerance)
            .unwrap_or_else(|e| panic!("volume: {e}")),
    )
}

fn assert_volume(rb: f64, x0: f64, len: f64, got: f64) {
    let want = analytic(rb, x0, len);
    let rel = (got - want).abs() / want;
    assert!(
        rel < 1e-4,
        "rb {rb} x0 {x0}: volume {got} vs analytic {want} (rel {rel:.2e})"
    );
}

#[test]
fn unequal_tees_serve_watertight_with_exact_volume() {
    for (rb, x0, len) in [
        (30.0, 0.0, 100.0),
        (30.0, -20.0, 120.0),
        (30.0, -100.0, 200.0),
        (10.0, 0.0, 100.0),
        (20.0, 0.0, 100.0),
        (35.0, 0.0, 100.0),
        (39.0, 0.0, 100.0),
    ] {
        let v = served_volume(rb, x0, len)
            .unwrap_or_else(|| panic!("rb {rb} x0 {x0}: the tee must serve"));
        assert_volume(rb, x0, len, v);
    }
}

#[test]
fn equal_radius_cross_serves_the_steinmetz_volume() {
    let v = served_volume(R, -100.0, 200.0).expect("the equal-radius cross must serve");
    // Two r = 40 cylinders through each other overlap in 16 r³ / 3.
    let want = std::f64::consts::PI * R * R * (H + 200.0) - 16.0 * R * R * R / 3.0;
    let rel = (v - want).abs() / want;
    assert!(rel < 1e-4, "volume {v} vs {want} (rel {rel:.2e})");
}

#[test]
fn equal_radius_tee_never_serves_its_operands_unmerged() {
    // The branch starting on the main axis: the PK-727 wrong serve.
    // It may decline until PK-732 builds its pinch vertices, but a
    // serve must carry the true volume, never the operands' sum.
    if let Some(v) = served_volume(R, 0.0, 100.0) {
        let sum = std::f64::consts::PI * R * R * (H + 100.0);
        assert!((v - sum).abs() > 1.0, "served the operands unmerged");
        // Overlap: ∫∫ √(R² − y²) over the branch disk = 2 ∫ (R² − y²) dy.
        let want = sum - 2.0 * (2.0 * R * R * R - 2.0 * R * R * R / 3.0);
        let rel = (v - want).abs() / want;
        assert!(rel < 1e-4, "volume {v} vs {want}");
    }
}
