//! A lone NURBS arc rim between two end faces (PK-769): the D profile's
//! semicircular top edge, filleted alone. Both ends are sharp corners on
//! the flat side, and the per-edge pipeline serves it (the chain surgery
//! leaves lone end-face edges per-edge). The strip ends on meridian
//! sections, so the removed volume is a half revolution of the fillet
//! section by Pappus.

use popkern_rhai::ScriptHost;

fn d_block(r: f64) -> String {
    format!(
        r#"
let m = model();
let b = m.profile(xy()).path(-20.0, 0.0).line_to(20.0, 0.0).arc_to([0.0, 20.0], [-20.0, 0.0]).close().sketch("d").extrude(blind(10.0));
b = b.fillet(edges().at_height(10.0).nurbs(), {r:?}, #{{ label: "rim" }});
m.assert_watertight();
m
"#
    )
}

#[test]
fn a_semicircular_rim_fillets_alone_and_runs_out_on_the_flat_side() {
    let (radius, h) = (20.0, 10.0);
    let pi = std::f64::consts::PI;
    for r in [1.0, 3.0] {
        let mut host = ScriptHost::new();
        let outcome = host
            .run(&d_block(r), "d.rhai")
            .unwrap_or_else(|e| panic!("r = {r}: {e}\n{}", e.payload));
        let mut model = outcome.model.expect("model");
        let part = model.to_part().expect("part");
        let brep = part.brep();
        let got = popkern_measure::body_volume(brep, &brep.tolerance).expect("volume");
        let a = (1.0 - pi / 4.0) * r * r;
        let x = r * (10.0 - 3.0 * pi) / (3.0 * (4.0 - pi));
        let want = 0.5 * pi * radius * radius * h - a * pi * (radius - x);
        let rel = (got - want).abs() / want;
        assert!(
            rel < 1e-4,
            "r = {r}: volume {got} vs {want} (rel {rel:.2e})"
        );
        occ_valid(&format!("d_{r}"), |step| {
            model.to_step(step).expect("STEP export");
        });
    }
}

/// OpenCASCADE reads the served body as one valid solid (skipped when
/// the OCP venv is absent).
fn occ_valid(name: &str, export: impl FnOnce(&std::path::Path)) {
    let home = std::path::PathBuf::from(std::env::var("HOME").unwrap_or_default());
    let python = home.join("venvs/ocp/bin/python");
    let tools = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tools/occ-oracle");
    if !python.exists() || !tools.exists() {
        eprintln!("[pk769] OCC oracle absent; validity not cross-checked");
        return;
    }
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("pk769");
    std::fs::create_dir_all(&dir).expect("dir");
    let step = dir.join(format!("{name}.step"));
    export(&step);
    let run = |script: &str| -> serde_json::Value {
        let out = std::process::Command::new(&python)
            .arg(tools.join(script))
            .arg(&step)
            .output()
            .expect("oracle");
        serde_json::from_slice(&out.stdout).expect("oracle json")
    };
    let check = run("occ_check.py");
    assert_eq!(
        check["passed"], true,
        "{name}: OCC rejects the STEP: {check}"
    );
    let solids = run("occ_solids.py");
    assert_eq!(
        solids["solids"].as_array().map(Vec::len),
        Some(1),
        "{name}: one solid: {solids}"
    );
}
