//! PK-527: the saddle setback's served bodies pass the OpenCASCADE
//! oracle. Both notched-cube saddles (constant and variable radius)
//! are exported to STEP. OpenCASCADE's BRepCheck must hold with zero
//! open free bounds.
//!
//! The `--self-interference` leg (BRepAlgoAPI_Check) takes about 10
//! minutes per body, so it is not a suite test. It measured both bodies
//! valid with no faults on 2026-09-25. To rerun it on the exported
//! files:
//! `~/venvs/ocp/bin/python tools/occ-oracle/occ_check.py
//! --self-interference target/tmp/saddle_setback_occ/*.step`
//!
//! The leg SKIPS when the OCP venv is absent, so CI without the oracle
//! stays green (the established oracle-battery pattern).

use popkern_blend::types::{FilletOptions, FilletRadius, FilletSpec};
use popkern_core::geometry::Point3;
use popkern_core::geometry::transform::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// A unit cube with its (+,+,+) corner notched out, and the saddle
/// vertex at `corner`.
fn notched_cube(corner: Point3) -> (BRep, popkern_core::topology::VertexHandle) {
    let tol = Tolerance::default();
    let notch = Transform3::translation(0.6, 0.6, 0.6).compose(&Transform3::scale_uniform(0.6));
    let tool = BRep::unit_cube(tol)
        .expect("cube")
        .transform(&notch)
        .expect("notch transform");
    let body = popkern_bool::subtraction(BRep::unit_cube(tol).expect("cube"), tool, tol, None)
        .expect("notched cube builds");
    let vertex = body
        .vertices
        .iter()
        .find(|(_, v)| v.position.distance_to(&corner) < 1e-6)
        .map(|(h, _)| h)
        .expect("saddle vertex");
    (body, vertex)
}

fn fillet_saddle(corner: Point3, radii: [FilletRadius; 3]) -> BRep {
    let (body, vertex) = notched_cube(corner);
    let mut edges: Vec<_> = body
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect();
    edges.sort();
    edges.dedup();
    let specs: Vec<FilletSpec> = edges
        .into_iter()
        .zip(radii)
        .map(|(edge, radius)| FilletSpec { edge, radius })
        .collect();
    let result = popkern_blend::fillet_edges(
        &body,
        &specs,
        &FilletOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect("saddle fillet serves");
    assert_eq!(
        result
            .corner_patches
            .iter()
            .filter(|p| p.vertex == vertex && p.setback_cap)
            .count(),
        16,
        "the saddle serves through the saddle setback (4 sectors x 4 pieces)"
    );
    result.brep
}

fn export_step_solid(brep: &BRep, path: &std::path::Path) -> Result<(), String> {
    let name = path
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "saddle".into());
    let doc = popkern_io::StepDocument {
        products: vec![popkern_io::StepProductNode {
            name: name.clone(),
            id: name.to_uppercase(),
            children: vec![],
            shapes: vec![std::sync::Arc::new(popkern_io::ShapeRep {
                kind: popkern_io::ShapeRepKind::AdvancedBrep,
                solids: vec![brep.clone()],
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
        .map_err(|e| format!("{e:?}"))
}

/// The oracle's verdict, or `None` when the OCP venv is absent.
fn occ_oracle_passes(path: &std::path::Path) -> Option<(bool, String)> {
    let python = std::path::Path::new(&std::env::var("HOME").unwrap_or_default())
        .join("venvs/ocp/bin/python");
    let script = std::path::Path::new("../../tools/occ-oracle/occ_check.py");
    if !python.exists() || !script.exists() {
        return None;
    }
    let output = std::process::Command::new(&python)
        .arg(script)
        .arg(path)
        .output()
        .ok()?;
    let report = String::from_utf8_lossy(&output.stdout).into_owned();
    Some((output.status.success(), report))
}

#[test]
fn saddle_setback_bodies_pass_the_occ_oracle() {
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("saddle_setback_occ");
    std::fs::create_dir_all(&dir).expect("oracle dir");
    let constant = [0.06, 0.08, 0.10].map(FilletRadius::Constant);
    let variable = [(0.05, 0.09), (0.07, 0.11), (0.08, 0.06)]
        .map(|(start, end)| FilletRadius::Variable { start, end });
    for (name, radii) in [("constant", constant), ("variable", variable)] {
        let brep = fillet_saddle(Point3::new(1.0, 0.6, 0.6), radii);
        let path = dir.join(format!("saddle_{name}.step"));
        export_step_solid(&brep, &path).expect("STEP export");
        let Some((passed, report)) = occ_oracle_passes(&path) else {
            eprintln!("[saddle-occ] OCP venv absent: skipping the oracle leg");
            return;
        };
        assert!(passed, "{name} saddle fails the OCC oracle: {report}");
    }
}
