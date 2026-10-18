//! The equal-radius tee's saddle fillet (PK-764, six-axis-arm §9.1 R1,
//! `docs/architecture/fillet-vanishing-ends.md`): an r 40 main cylinder
//! and an r 40 branch from its axis meet in two half-ellipse saddle
//! chains S− → seam crossing → S+. At the singular points S± the
//! cylinders are tangent, so the fillet strips run out to zero width.
//! Each chain is built whole: its seam-crossing junction is exact and
//! its strips end on the tangent points, with no corner blend there.
//!
//! - It serves watertight at r = 2, 4 and 8 with no failure. At
//!   r = 2 and 4 every seam certifies G1 (r = 8: PK-767).
//! - OpenCASCADE reads one valid solid.
//! - The concave fillet adds material along a fixed saddle, so
//!   (V − V₀) / r² agrees across the radii. It approaches the
//!   first-order section integral as r shrinks: 21.63, the wedge section
//!   r²(cot(β/2) − (π − β)/2) integrated along both chains. OCC's
//!   adaptive volumes are used: the added volume is 1e-4 of the body,
//!   below popkern's mesh-volume accuracy and OCC's default rule.

use popkern_rhai::ScriptHost;

fn tee(fillet: Option<f64>) -> String {
    let line = fillet
        .map(|r| {
            format!(
                r#"b = b.fillet(edges().between("cylinder", "cylinder"), {r:?}, #{{ label: "saddle" }});"#
            )
        })
        .unwrap_or_default();
    format!(
        r#"
let m = model();
let main = m.profile(xy()).circle(40.0).sketch("main").extrude(blind(120.0));
let branch = m.tool(|t| t.profile(yz_at(0.0)).circle(40.0).at(0.0, 60.0).sketch("branch").extrude(blind(100.0)));
let b = main.union(branch, #{{ label: "tee" }});
{line}
m.assert_watertight();
m
"#
    )
}

/// OCC's exact volume of the served body, after OCC's validity check;
/// `None` when the OCP venv is absent.
fn occ_volume(source: &str, name: &str) -> Option<f64> {
    let home = std::path::PathBuf::from(std::env::var("HOME").unwrap_or_default());
    let python = home.join("venvs/ocp/bin/python");
    let tools = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tools/occ-oracle");
    if !python.exists() || !tools.exists() {
        eprintln!("[pk764] OCC oracle absent; volumes not cross-checked");
        return None;
    }
    let mut host = ScriptHost::new();
    let mut model = host
        .run(source, "tee.rhai")
        .unwrap_or_else(|e| panic!("{}", e.render()))
        .model
        .expect("model");
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("pk764");
    std::fs::create_dir_all(&dir).expect("dir");
    let step = dir.join(format!("{name}.step"));
    model.to_step(&step).expect("STEP export");
    let run = |script: &str, extra: &[&str]| -> serde_json::Value {
        let out = std::process::Command::new(&python)
            .arg(tools.join(script))
            .arg(&step)
            .args(extra)
            .output()
            .expect("oracle");
        serde_json::from_slice(&out.stdout).expect("oracle json")
    };
    let check = run("occ_check.py", &[]);
    assert_eq!(
        check["passed"], true,
        "{name}: OCC rejects the STEP: {check}"
    );
    let solids = run("occ_solids.py", &["--adaptive"]);
    let solids = solids["solids"].as_array().expect("solids");
    assert_eq!(solids.len(), 1, "{name}: one solid");
    solids[0]["volume"].as_f64()
}

#[test]
fn the_equal_radius_saddle_fillets_with_strips_that_run_out_at_the_singular_points() {
    let mut added = Vec::new();
    let v0 = occ_volume(&tee(None), "tee");
    for r in [2.0, 4.0, 8.0] {
        let source = tee(Some(r));
        let mut host = ScriptHost::new();
        let outcome = host
            .run(&source, "tee.rhai")
            .unwrap_or_else(|e| panic!("r {r}: {}", e.render()));
        let model = outcome.model.expect("model");
        let report = popkern_rhai::report_value(&model);
        assert_eq!(report["failures"], 0, "r {r}: {report}");
        let saddle = report["features"]
            .as_array()
            .and_then(|f| f.iter().find(|row| row["label"] == "saddle"))
            .expect("saddle row");
        let achieved = &saddle["blend"]["continuity"]["achieved"];
        eprintln!("[pk764] r {r}: continuity achieved {achieved}");
        // At r = 2 and r = 4 every seam certifies G1 (r = 4 since the
        // 2026-10 station and seam work; re-measured 2026-10-18). At
        // r = 8 the blend-to-support seams part 9.6e-5 from their
        // supports near the vanishing tip, over the 1e-5 G0 gate, while
        // the walker's rail midpoints sit at 8.5e-6: PK-767.
        if r <= 4.0 {
            assert_eq!(achieved, "G1", "r {r}: {saddle}");
        }
        if let (Some(v0), Some(v)) = (v0, occ_volume(&source, &format!("fillet_{r}"))) {
            assert!(
                v > v0,
                "r {r}: a concave fillet adds material ({v} vs {v0})"
            );
            added.push((v - v0) / (r * r));
        }
    }
    if let (Some(lo), Some(hi)) = (
        added.iter().copied().reduce(f64::min),
        added.iter().copied().reduce(f64::max),
    ) {
        eprintln!("[pk764] (V - V0)/r^2 across r = 2, 4, 8: {added:?}");
        assert!(
            (hi - lo) / hi < 0.03,
            "section scaling inconsistent: {added:?}"
        );
        // Toward the first-order section integral as r shrinks.
        assert!(
            (added[0] - 21.633).abs() / 21.633 < 0.01,
            "r = 2: {added:?}"
        );
        assert!(
            added.windows(2).all(|w| w[1] < w[0]),
            "smooth in r: {added:?}"
        );
    }
}
