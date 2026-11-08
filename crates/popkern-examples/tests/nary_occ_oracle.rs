//! N-ary booleans on the partition against OpenCASCADE's own fold of
//! the same operands (PK-1048, design 10): union, common and the N-body
//! symmetric difference (exactly one body) of curved solids, compared
//! by volume. OCC has no exactly-one of its own; `occ_fuse.py
//! --op exactly-one` cuts each operand by the fuse of the others, whose
//! pieces are disjoint. The two kernels share no boolean code.
//!
//! The leg SKIPS when the OCP venv is absent (CI without the oracle
//! stays green).

use popkern_bool::{BoolOptions, NonManifoldPolicy, SetExpr, nary_boolean_bodies};
use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn export(solids: Vec<BRep>, path: &std::path::Path) {
    let name = path
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_default();
    let doc = popkern_io::StepDocument {
        products: vec![popkern_io::StepProductNode {
            name: name.clone(),
            id: name.to_uppercase(),
            children: vec![],
            shapes: vec![std::sync::Arc::new(popkern_io::ShapeRep {
                kind: popkern_io::ShapeRepKind::AdvancedBrep,
                solids,
                context: popkern_io::RepresentationContext::default(),
                non_brep_items: vec![],
            })],
            metadata: popkern_io::ProductMetadata::default(),
        }],
        roots: vec![0],
        units: popkern_io::UnitContext::default(),
        opaque_bag: popkern_io::OpaqueBag::default(),
    };
    popkern_io::export_step_file(&doc, path, &popkern_io::StepExportOptions::default())
        .expect("STEP export");
}

/// Run `script` with `args` in order; `None` without the OCP venv.
fn oracle(script: &str, args: &[&std::ffi::OsStr]) -> Option<serde_json::Value> {
    let python = std::path::Path::new(&std::env::var("HOME").unwrap_or_default())
        .join("venvs/ocp/bin/python");
    let tools = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tools/occ-oracle");
    if !python.exists() || !tools.exists() {
        eprintln!("[nary-oracle] OCC oracle absent; volumes not cross-checked");
        return None;
    }
    let out = std::process::Command::new(&python)
        .arg(tools.join(script))
        .args(args)
        .output()
        .expect("oracle");
    Some(serde_json::from_slice(&out.stdout).expect("oracle json"))
}

/// The served bodies of `expr` over `bodies` (separate manifold bodies
/// where the kept regions touch at edges or vertices) against OCC's `op`
/// fold of the same operands: their volumes agree to 1e-6 relative.
fn agrees(name: &str, bodies: &[BRep], expr: &SetExpr, op: &str) {
    let tol = Tolerance::default();
    let refs: Vec<&BRep> = bodies.iter().collect();
    let options = BoolOptions {
        non_manifold: NonManifoldPolicy::SplitToManifold,
    };
    let served =
        nary_boolean_bodies(&refs, expr, tol, &options).unwrap_or_else(|e| panic!("{name}: {e}"));
    for body in &served {
        assert!(body.validate().is_ok(), "{name}: {:?}", body.validate());
    }
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("nary_oracle");
    std::fs::create_dir_all(&dir).expect("dir");
    let step = dir.join(format!("{name}.step"));
    export(served, &step);
    let mut operands = Vec::new();
    for (k, b) in bodies.iter().enumerate() {
        let path = dir.join(format!("{name}_operand_{k}.step"));
        export(vec![b.clone()], &path);
        operands.push(path);
    }
    let Some(check) = oracle("occ_check.py", &[step.as_os_str()]) else {
        return;
    };
    assert_eq!(
        check["passed"], true,
        "{name}: OCC rejects the STEP: {check}"
    );
    let solids =
        oracle("occ_solids.py", &[step.as_os_str(), "--adaptive".as_ref()]).expect("oracle");
    let served_volume: f64 = solids["solids"]
        .as_array()
        .expect("solids")
        .iter()
        .map(|r| r["volume"].as_f64().expect("volume"))
        .sum();
    let mut args: Vec<&std::ffi::OsStr> = vec!["--op".as_ref(), op.as_ref()];
    args.extend(operands.iter().map(|p| p.as_os_str()));
    let folded = oracle("occ_fuse.py", &args).expect("oracle");
    assert_eq!(folded["done"], true, "{name}: OCC {op} {folded}");
    let occ = folded["volume"].as_f64().expect("volume");
    eprintln!("[nary-oracle] {name}: served {served_volume:.9} occ {occ:.9}");
    assert!(
        (served_volume - occ).abs() <= 1e-6 * occ.max(1.0),
        "{name}: served {served_volume} occ {op} {occ}"
    );
}

/// Rods of different radii along tilted axes through near the origin.
fn rods(n: usize) -> Vec<BRep> {
    let tol = Tolerance::default();
    (0..n)
        .map(|k| {
            BRep::cylinder(0.4 + 0.03 * k as f64, 3.0, tol)
                .and_then(|c| c.transform(&Transform3::translation(0.0, 0.017 * k as f64, -1.5)))
                .and_then(|c| c.transform(&Transform3::rotation_y(std::f64::consts::FRAC_PI_2)))
                .and_then(|c| c.transform(&Transform3::rotation_z(0.7 * k as f64 + 0.05)))
                .expect("rod")
        })
        .collect()
}

fn sphere_and_rod() -> Vec<BRep> {
    let tol = Tolerance::default();
    let sphere = BRep::sphere(1.0, tol).expect("sphere");
    // Across the sphere's axis: a rod holding both poles is PK-1055.
    let rod = BRep::cylinder(0.6, 3.0, tol)
        .and_then(|c| c.transform(&Transform3::translation(0.1, 0.05, -1.5)))
        .and_then(|c| c.transform(&Transform3::rotation_y(std::f64::consts::FRAC_PI_2)))
        .expect("rod");
    vec![sphere, rod]
}

#[test]
fn sphere_and_rod_symmetric_difference_matches_occ() {
    agrees(
        "sphere_rod_xor",
        &sphere_and_rod(),
        &SetExpr::ExactlyOne,
        "exactly-one",
    );
}

#[test]
fn three_rods_exactly_one_matches_occ() {
    agrees(
        "rods3_exactly_one",
        &rods(3),
        &SetExpr::ExactlyOne,
        "exactly-one",
    );
}

#[test]
fn four_rods_union_matches_occ() {
    agrees("rods4_union", &rods(4), &SetExpr::Any, "fuse");
}

#[test]
fn four_rods_common_matches_occ() {
    agrees("rods4_common", &rods(4), &SetExpr::All, "common");
}
