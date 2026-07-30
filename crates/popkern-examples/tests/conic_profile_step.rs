//! Conic blend profile STEP/OCC acceptance (PK-189, requirements
//! 2.2 acceptance 273): a conic-profile blend must serve
//! watertight, export to STEP, and survive the OCC oracle. The
//! oracle leg follows the corpus pattern — it skips with a notice
//! when the OCP venv is absent so CI without the oracle stays
//! green.

use popkern_part::prelude::*;

#[test]
fn conic_profile_plate_serves_exports_and_survives_occ() {
    let plate = Profile::on(Plane::XY)
        .rect(60.0, 40.0)
        .extrude(8.0)
        .unwrap();
    let part = plate
        .fillet_with(
            |recipe| recipe.radius(3.0).conic(0.35),
            |e| {
                e.straight()
                    .direction(Dir::X)
                    .chord_length_between(55.0, 65.0)
                    .at_height(8.0)
            },
        )
        .unwrap();
    let brep = part.brep();
    assert!(
        brep.validate_all().is_empty(),
        "conic plate must validate: {:?}",
        brep.validate_all()
    );
    let opens = part.open_edge_count().expect("open edge count");
    assert_eq!(opens, 0, "conic plate must serve watertight");

    // The conic (rho = 0.35, flatter than circular √2−1) removes
    // MORE material than the circular fillet; sanity-pin the volume
    // between the circular-fillet value and the full chamfer wedge.
    let volume = popkern_measure::body_volume(brep, &popkern_core::tolerance::Tolerance::default())
        .expect("volume");
    let base = 60.0 * 40.0 * 8.0;
    let circular_removal = 2.0 * (1.0 - std::f64::consts::PI / 4.0) * 3.0 * 3.0 * 60.0;
    let chamfer_removal = 2.0 * 0.5 * 3.0 * 3.0 * 60.0;
    assert!(
        volume < base - circular_removal && volume > base - chamfer_removal,
        "conic volume {volume} must sit between the circular fillet \
         ({}) and the chamfer wedge ({})",
        base - circular_removal,
        base - chamfer_removal
    );

    let out_dir = std::path::Path::new("../../examples/output/fixtures");
    std::fs::create_dir_all(out_dir).expect("fixture dir");
    let path = out_dir.join("conic_profile_plate.step");
    part.export_step(&path).expect("conic STEP export");

    let oracle_python = std::path::Path::new(&std::env::var("HOME").unwrap_or_default())
        .join("venvs/ocp/bin/python");
    let oracle_script = std::path::Path::new("../../tools/occ-oracle/occ_check.py");
    if !oracle_python.exists() || !oracle_script.exists() {
        eprintln!("[conic-occ] OCP oracle unavailable — leg skipped");
        return;
    }
    let output = std::process::Command::new(&oracle_python)
        .arg(oracle_script)
        .arg(&path)
        .output()
        .expect("oracle invocation");
    let report = String::from_utf8_lossy(&output.stdout);
    assert!(
        report.contains("\"passed\": true"),
        "OCC rejected the conic-profile STEP: {}",
        report.chars().take(400).collect::<String>()
    );
}
