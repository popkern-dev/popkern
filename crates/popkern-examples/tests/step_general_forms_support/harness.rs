//! The round trip `step_general_forms` runs: a body written as a
//! one-product document, read back under an option, and compared by
//! census, solid volumes and OpenCASCADE's reading.

use std::sync::Arc;

use popkern_core::topology::{BRep, RegionKind};
use popkern_io::step::diagnostics::Severity;
use popkern_io::{
    OpaqueBag, ProductMetadata, RepresentationContext, ShapeRep, ShapeRepKind, StepDocument,
    StepImportOptions, StepNonManifold, StepProductNode, UnitContext, import_step,
};

pub(crate) fn document(brep: &BRep, name: &str) -> StepDocument {
    StepDocument {
        products: vec![StepProductNode {
            name: name.into(),
            id: name.into(),
            children: Vec::new(),
            shapes: vec![Arc::new(ShapeRep {
                kind: ShapeRepKind::AdvancedBrep,
                solids: vec![brep.clone()],
                context: RepresentationContext::default(),
                non_brep_items: Vec::new(),
            })],
            metadata: ProductMetadata::default(),
        }],
        roots: vec![0],
        units: UnitContext::default(),
        opaque_bag: OpaqueBag::default(),
    }
}

/// Every body the file holds under `non_manifold`, with no error
/// diagnostics.
pub(crate) fn bodies(data: &[u8], non_manifold: StepNonManifold, name: &str) -> Vec<BRep> {
    let imported = import_step(
        data,
        &StepImportOptions {
            non_manifold,
            ..Default::default()
        },
    )
    .expect("import");
    let errors: Vec<_> = imported
        .diagnostics
        .entries()
        .iter()
        .filter(|e| e.severity == Severity::Error)
        .collect();
    assert!(errors.is_empty(), "{name}: {errors:?}");
    imported
        .document
        .products
        .iter()
        .flat_map(|p| &p.shapes)
        .flat_map(|s| s.solids.iter().cloned())
        .collect()
}

/// What the test compares: entity counts, then for each region kind
/// (outside included) how many shells hold faces, wires and acorns.
/// Vertex, edge and face counts, then one row per shell: its region
/// kind and how many faces, wire edges and acorns it holds.
pub(crate) type Census = ((usize, usize, usize), Vec<(String, usize, usize, usize)>);

pub(crate) fn census(body: &BRep) -> Census {
    let counts = (body.vertices.len(), body.edges.len(), body.faces.len());
    let mut rows: Vec<(String, usize, usize, usize)> = body
        .shells
        .iter()
        .map(|(h, s)| {
            let kind = if s.region == body.outside {
                "outside".to_string()
            } else {
                format!("{:?}", body.region_kind(s.region))
            };
            let faces = body
                .faces
                .values()
                .filter(|f| f.back == h || f.front == h)
                .count();
            (
                kind,
                faces,
                s.wire_edges.len(),
                usize::from(s.acorn.is_some()),
            )
        })
        .collect();
    rows.sort();
    (counts, rows)
}

pub(crate) fn solid_volumes(body: &BRep) -> Vec<f64> {
    let mut v: Vec<f64> = popkern_measure::region_volumes(body, &body.tolerance)
        .expect("volumes")
        .into_iter()
        .filter(|(r, _)| body.region_kind(*r) == Some(RegionKind::Solid))
        .map(|(_, v)| v)
        .collect();
    v.sort_by(f64::total_cmp);
    v
}

/// OpenCASCADE's solid volumes and distinct face, edge and vertex
/// counts from `data`; `None` without the OCP venv.
pub(crate) fn occ_read(data: &[u8], name: &str) -> Option<(Vec<f64>, [u64; 3])> {
    let home = std::path::PathBuf::from(std::env::var("HOME").unwrap_or_default());
    let python = home.join("venvs/ocp/bin/python");
    let oracle = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tools/occ-oracle/occ_solids.py");
    if !python.exists() || !oracle.exists() {
        return None;
    }
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("step_general_forms");
    std::fs::create_dir_all(&dir).expect("dir");
    let path = dir.join(format!("{name}.step"));
    std::fs::write(&path, data).expect("write");
    let out = std::process::Command::new(&python)
        .arg(&oracle)
        .arg(&path)
        .output()
        .expect("oracle");
    let report: serde_json::Value = serde_json::from_slice(&out.stdout).expect("json");
    let mut volumes: Vec<f64> = report["solids"]
        .as_array()
        .expect("solids")
        .iter()
        .map(|s| s["volume"].as_f64().expect("volume"))
        .collect();
    volumes.sort_by(f64::total_cmp);
    let count = |k: &str| report["counts"][k].as_u64().expect("count");
    Some((volumes, [count("faces"), count("edges"), count("vertices")]))
}

pub(crate) fn assert_close(name: &str, ours: &[f64], want: &[f64]) {
    assert_eq!(ours.len(), want.len(), "{name}: {ours:?} vs {want:?}");
    for (o, w) in ours.iter().zip(want) {
        assert!(
            (o - w).abs() < 1e-6 * w.abs().max(1.0),
            "{name}: {o} vs {w}"
        );
    }
}
