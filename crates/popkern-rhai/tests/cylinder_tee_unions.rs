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
//!   at the crossings.
//! - PK-732: with the branch end on the main axis, its end circle
//!   TOUCHES the main cylinder at the two singular points; the touch
//!   arrived as two edge-face roots 6e-4 apart (the tolerance band
//!   around a tangency), the circle split twice there and the walk
//!   stranded the records. One contact is one root: union, subtract
//!   and intersect serve for the axis tee and the x0 = −20 tee.
//!
//! Volumes are checked against the analytic value: the overlap of a
//! branch of radius `rb` through the main cylinder is
//! `∫ 2√(R²−y²) · 2√(rb²−y²) dy` over `|y| < rb`.

use popkern_rhai::ScriptHost;

const R: f64 = 40.0;
const H: f64 = 120.0;

fn tee(rb: f64, x0: f64, len: f64) -> String {
    tee_op("union", rb, x0, len)
}

/// The tee under boolean `op` (`union`, `subtract` or `intersect`).
fn tee_op(op: &str, rb: f64, x0: f64, len: f64) -> String {
    format!(
        r#"
let m = model();
let main = m.profile(xy()).circle({R:?}).sketch("main").extrude(blind({H:?}));
let branch = m.tool(|t| t.profile(yz_at({x0:?})).circle({rb:?}).at(0.0, 60.0).sketch("branch").extrude(blind({len:?})));
let b = main.{op}(branch, #{{ label: "tee" }});
m.assert_watertight();
m
"#
    )
}

/// The volume of main ∩ branch for a branch spanning `x ∈ [x0, x0 + len]`
/// (`x0 + len ≥ R`): at height offset `y` the main cylinder's chord is
/// `|x| ≤ hc`, of which the branch covers `x ≥ x0`.
fn overlap(rb: f64, x0: f64) -> f64 {
    let n = 200_000;
    let mut sum = 0.0;
    for i in 0..n {
        let y = -rb + 2.0 * rb * (f64::from(i) + 0.5) / f64::from(n);
        let hc = (R * R - y * y).max(0.0).sqrt();
        let depth = (hc - x0.max(-hc)).max(0.0);
        sum += depth * 2.0 * (rb * rb - y * y).max(0.0).sqrt() * 2.0 * rb / f64::from(n);
    }
    sum
}

/// The volume of main ∪ branch.
fn analytic(rb: f64, x0: f64, len: f64) -> f64 {
    std::f64::consts::PI * (R * R * H + rb * rb * len) - overlap(rb, x0)
}

/// Run the tee: `Some(volume)` when it serves (one watertight body),
/// `None` when the kernel declines.
fn served_volume(rb: f64, x0: f64, len: f64) -> Option<f64> {
    served_volume_of(&tee(rb, x0, len))
}

/// The served body's volume by OpenCASCADE's exact integration of its
/// STEP export (`None` when the OCP venv is absent): popkern's
/// `body_volume` integrates a fine mesh, accurate to about 1e-5 on
/// these tangent-cornered trims.
fn occ_volume_of(script: &str, name: &str) -> Option<f64> {
    let home = std::path::PathBuf::from(std::env::var("HOME").unwrap_or_default());
    let python = home.join("venvs/ocp/bin/python");
    let oracle = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tools/occ-oracle/occ_solids.py");
    if !python.exists() || !oracle.exists() {
        return None;
    }
    let mut host = ScriptHost::new();
    let mut model = host.run(script, "tee.rhai").ok()?.model?;
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("pk732");
    std::fs::create_dir_all(&dir).expect("dir");
    let step = dir.join(format!("{name}.step"));
    model.to_step(&step).expect("STEP export");
    let out = std::process::Command::new(&python)
        .arg(&oracle)
        .arg(&step)
        .output()
        .expect("oracle");
    let report: serde_json::Value = serde_json::from_slice(&out.stdout).expect("oracle json");
    let solids = report["solids"].as_array().expect("solids");
    assert_eq!(solids.len(), 1, "{name}: OCC reads one solid");
    // OCC's own validity check: one shell, no free boundaries.
    let check = std::process::Command::new(&python)
        .arg(oracle.with_file_name("occ_check.py"))
        .arg(&step)
        .output()
        .expect("occ_check");
    let verdict: serde_json::Value = serde_json::from_slice(&check.stdout).expect("check json");
    assert_eq!(
        verdict["passed"], true,
        "{name}: OCC rejects the STEP: {verdict}"
    );
    solids[0]["volume"].as_f64()
}

fn served_volume_of(script: &str) -> Option<f64> {
    let mut host = ScriptHost::new();
    let outcome = host.run(script, "tee.rhai").ok()?;
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
fn equal_radius_tees_serve_every_boolean_with_exact_volume() {
    // The branch end on the main axis (its end circle touches the main
    // cylinder at the singular points; PK-727 served the operands
    // unmerged) and inside it at x0 = −20 (the crossings lie on the
    // branch wall).
    let main = std::f64::consts::PI * R * R * H;
    for x0 in [0.0, -20.0] {
        let len = 100.0 - x0;
        let branch = std::f64::consts::PI * R * R * len;
        let both = overlap(R, x0);
        for (op, want) in [
            ("union", main + branch - both),
            ("subtract", main - both),
            ("intersect", both),
        ] {
            let script = tee_op(op, R, x0, len);
            let got = served_volume_of(&script)
                .unwrap_or_else(|| panic!("{op} x0 {x0}: the equal-radius tee must serve"));
            let rel = (got - want).abs() / want;
            assert!(
                rel < 1e-4,
                "{op} x0 {x0}: mesh volume {got} vs {want} (rel {rel:.2e})"
            );
            if let Some(exact) = occ_volume_of(&script, &format!("{op}_{}", -x0)) {
                let rel = (exact - want).abs() / want;
                assert!(
                    rel < 1e-6,
                    "{op} x0 {x0}: OCC volume {exact} vs {want} (rel {rel:.2e})"
                );
            }
        }
    }
    // The axis tee's overlap in closed form: ∫∫ √(R² − y²) over the
    // branch disk = 2 ∫ (R² − y²) dy = 8R³/3.
    let closed = 8.0 * R * R * R / 3.0;
    assert!((overlap(R, 0.0) - closed).abs() / closed < 1e-6);
}

/// A near-tangent tee's intersection edges keep their pcurves on the
/// edge (PK-885): the projection fit's refinement quit at 44 spans
/// when a doubling improved the steep end spans by less than half,
/// serving the r 39 branch's edges 1.03e-4 off — over the model's
/// max_imprecision, recorded at the ceiling until honest recording
/// declined the union.
#[test]
fn near_tangent_tee_edges_keep_their_pcurves_on_the_edge() {
    let mut host = ScriptHost::new();
    let mut model = host
        .run(&tee(39.0, 0.0, 100.0), "tee.rhai")
        .expect("tee serves")
        .model
        .expect("model");
    let body = model.to_brep().expect("body");
    let worst = body
        .half_edges
        .keys()
        .filter_map(|he| body.measure_half_edge_pcurve_deviation(he))
        .map(|m| m.max_deviation)
        .fold(0.0_f64, f64::max);
    assert!(worst <= 1e-5, "worst pcurve image {worst:.3e}");
}
