//! Interchange conditioning of setback vertex caps (PK-524).
//!
//! The dual cross-jet corner cells of a setback cap are exact
//! rational bi-11 patches whose weights span 1e-14 .. 2 (their
//! squared quadrant weights vanish at the cell corner). Written as-is,
//! OpenCASCADE's STEP import inflated the triangular-pyramid cap's
//! rim-quarter edge tolerances to 5.5e-2 and its boolean pre-check
//! flagged neighbouring corners as self-interfering, while every
//! curve, pcurve and vertex was exact in popkern. The exporter now
//! writes such a patch as a measured non-rational fit over the same
//! chart.
//!
//! Gates: no B-spline surface in the written file carries a weight
//! spread below 1e-8, the file round-trips through popkern's importer
//! without errors, and — when the local OCC oracle venv exists — the
//! oracle passes WITH its self-interference check.

use popkern_blend::{ContinuityLevel, FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_io::{
    OpaqueBag, ProductMetadata, RepresentationContext, ShapeRep, ShapeRepKind, StepDocument,
    StepExportOptions, StepImportOptions, StepProductNode, UnitContext, export_step, import_step,
};
use std::path::Path;
use std::sync::Arc;

fn apex_cap(n: usize) -> BRep {
    let tol = Tolerance::default();
    let pyramid = BRep::n_gonal_pyramid(n, 1.0, 1.2, tol).expect("pyramid");
    let apex = pyramid
        .vertices
        .iter()
        .max_by(|a, b| a.1.position.z.total_cmp(&b.1.position.z))
        .map(|(h, _)| h)
        .expect("apex");
    let specs: Vec<FilletSpec> = pyramid
        .edges
        .iter()
        .filter_map(|(eh, edge)| {
            let he = &pyramid.half_edges[edge.half_edge];
            (he.origin == apex || pyramid.half_edges[he.twin].origin == apex).then_some(eh)
        })
        .map(|edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.16),
        })
        .collect();
    let options = FilletOptions {
        continuity: ContinuityLevel::G2,
        ..Default::default()
    };
    fillet_edges(&pyramid, &specs, &options, tol, None)
        .expect("apex fillet serves")
        .brep
}

fn step_bytes(brep: BRep, name: &str) -> Vec<u8> {
    let document = StepDocument {
        products: vec![StepProductNode {
            name: name.into(),
            id: name.to_uppercase(),
            children: vec![],
            shapes: vec![Arc::new(ShapeRep {
                kind: ShapeRepKind::AdvancedBrep,
                solids: vec![brep],
                context: RepresentationContext::default(),
                non_brep_items: vec![],
            })],
            metadata: ProductMetadata::default(),
        }],
        roots: vec![0],
        units: UnitContext::default(),
        opaque_bag: OpaqueBag::default(),
    };
    export_step(&document, &StepExportOptions::default()).expect("export")
}

fn oracle_with_self_interference(step: &Path) -> Option<bool> {
    let python = Path::new(&std::env::var("HOME").unwrap_or_default()).join("venvs/ocp/bin/python");
    let check = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tools/occ-oracle/occ_check.py");
    if !(python.exists() && check.exists()) {
        return None;
    }
    let output = std::process::Command::new(&python)
        .arg(&check)
        .arg("--self-interference")
        .arg(step)
        .output()
        .expect("oracle");
    Some(String::from_utf8_lossy(&output.stdout).contains("\"passed\": true"))
}

#[test]
fn setback_caps_write_well_conditioned_surfaces_that_occ_accepts() {
    for n in [3usize, 4, 5] {
        let name = format!("apex_cap_n{n}");
        let bytes = step_bytes(apex_cap(n), &name);
        let imported = import_step(&bytes, &StepImportOptions::default()).expect("import");
        assert!(
            imported.diagnostics.errors().next().is_none(),
            "n={n}: round-trip import errors: {:?}",
            imported.diagnostics.entries()
        );
        let body = &imported.document.products[0].shapes[0].solids[0];
        for (face, data) in body.faces.iter() {
            let Some(SurfaceGeometry::Nurbs(nurbs)) = data.surface.as_ref().map(|s| &s.geometry)
            else {
                continue;
            };
            let (lo, hi) = nurbs
                .weights
                .iter()
                .fold((f64::INFINITY, 0.0_f64), |(a, b), &w| (a.min(w), b.max(w)));
            assert!(
                lo / hi >= 1.0e-8,
                "n={n}: {face:?} was written with weight spread {:.1e}..{:.1e}",
                lo,
                hi
            );
        }
        let path = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("{name}.step"));
        std::fs::write(&path, &bytes).expect("write");
        if let Some(passed) = oracle_with_self_interference(&path) {
            assert!(passed, "n={n}: OCC oracle with self-interference must pass");
        }
    }
}
