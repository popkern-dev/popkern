//! §5.3 item #2 — STEP round-trip regression on a Gregory-carrying BRep.
//!
//! Closes the remaining deliverable on item #2: "STEP export regression:
//! filleted mixed-radius cube → write_step → read_step → validate
//! invariants preserved." STEP AP203/242 has no Gregory entity, so
//! export goes through `SurfaceGeometry::to_nurbs()` (which for a
//! Gregory `[0,1]²` range calls `to_nurbs_exact` — the exact rational
//! rectangular-Bezier form at bidegree (7,7) for cubic or (9,9) for
//! quintic). After round-trip, Gregory faces return as NURBS faces;
//! this test verifies solid/face counts and sampled geometric
//! fidelity are preserved.

use std::sync::Arc;

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{
    BlendResult, CornerMethod, CornerPatchInfo, FilletOptions, FilletRadius, FilletSpec,
    fillet_edges,
};
use popkern_io::{
    OpaqueBag, ProductMetadata, RepresentationContext, ShapeRep, ShapeRepKind, StepDocument,
    StepExportOptions, StepImportOptions, StepProductNode, UnitContext, export_step, import_step,
};

fn mixed_radius_gregory_cube_result() -> BlendResult {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let incident_edges: Vec<_> = cube
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect();
    let radii = [0.08, 0.10, 0.12];
    let specs: Vec<FilletSpec> = incident_edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let opts = FilletOptions {
        use_rational_gregory: true,
        ..FilletOptions::default()
    };
    fillet_edges(&cube, &specs, &opts, Tolerance::default(), None)
        .expect("mixed-radius Gregory cube must build")
}

fn mixed_radius_gregory_cube() -> BRep {
    mixed_radius_gregory_cube_result().brep
}

struct CornerPatchFaceSamples {
    samples: Vec<(Point3, f64, f64)>,
}

fn sample_gregory_method_corner_patch_faces(
    brep: &BRep,
    corner_patches: &[CornerPatchInfo],
) -> Vec<CornerPatchFaceSamples> {
    let mut faces = Vec::new();
    for info in corner_patches {
        if info.method != CornerMethod::Gregory {
            continue;
        }
        let Some(sof) = brep
            .faces
            .get(info.face)
            .and_then(|face| face.surface.as_ref())
        else {
            continue;
        };
        let (SurfaceGeometry::Gregory(_) | SurfaceGeometry::Nurbs(_)) = &sof.geometry else {
            continue;
        };
        let mut samples = Vec::new();
        // Sample on a 3×3 interior UV grid — avoid corners where
        // ε-regularization introduces ~1e-14 noise.
        for i in 1..=3 {
            for j in 1..=3 {
                let u = i as f64 / 4.0;
                let v = j as f64 / 4.0;
                samples.push((sof.geometry.point_at(u, v), u, v));
            }
        }
        faces.push(CornerPatchFaceSamples { samples });
    }
    faces
}

/// STEP round-trip on a mixed-radius Gregory fillet. Post-round-trip
/// the BRep must:
/// (1) re-import without errors,
/// (2) carry ≥1 solid with ≥ the original face count,
/// (3) validate without structural errors.
#[test]
fn gregory_fillet_step_roundtrip_preserves_topology() {
    let brep = mixed_radius_gregory_cube();
    let orig_face_count = brep.faces.len();

    let doc = StepDocument {
        products: vec![StepProductNode {
            name: "GregoryFilletedCube".into(),
            id: "GFC001".into(),
            children: vec![],
            shapes: vec![Arc::new(ShapeRep {
                kind: ShapeRepKind::AdvancedBrep,
                solids: vec![brep.clone()],
                context: RepresentationContext::default(),
                non_brep_items: vec![],
            })],
            metadata: ProductMetadata::default(),
        }],
        roots: vec![0],
        units: UnitContext::default(),
        opaque_bag: OpaqueBag::default(),
    };
    let bytes = export_step(&doc, &StepExportOptions::default())
        .expect("export_step must succeed on Gregory-carrying BRep");
    assert!(
        bytes.len() > 1000,
        "Gregory-faced BRep STEP encoding should produce non-trivial output; got {} bytes",
        bytes.len()
    );

    let import_result =
        import_step(&bytes, &StepImportOptions::default()).expect("import_step must succeed");
    let reimported = &import_result.document.products[0].shapes[0].solids[0];

    assert!(
        reimported.faces.len() >= orig_face_count,
        "re-imported BRep must have ≥ {orig_face_count} faces; got {}",
        reimported.faces.len()
    );
    assert_eq!(
        reimported.solids.len(),
        1,
        "re-imported BRep must have 1 solid"
    );

    // Structural validation: re-imported BRep must not have broken
    // cycles, broken twins, or non-manifold edges.
    use popkern_core::error::ValidationError;
    let errors = reimported.validate_all();
    let structural: Vec<_> = errors
        .iter()
        .filter(|e| {
            matches!(
                e,
                ValidationError::BrokenCycle { .. }
                    | ValidationError::BrokenTwin { .. }
                    | ValidationError::NonManifoldEdge { .. }
            )
        })
        .collect();
    assert!(
        structural.is_empty(),
        "re-imported BRep has structural errors: {structural:?}"
    );
}

/// The Gregory-method cap faces in the original BRep are converted to
/// STEP NURBS geometry on export. Geometric fidelity must be preserved:
/// sample the original cap faces at interior UV points, locate the
/// corresponding face in the re-imported BRep, and assert the position
/// deviation is within tolerance.
///
/// Older Gregory wrappers use `to_nurbs_exact`; newer macro-cap rungs
/// may already be NURBS before export. In both cases STEP should carry
/// the same-UV surface samples to < 1e-6 after parser re-encoding.
#[test]
fn gregory_face_geometry_survives_step_roundtrip() {
    let result = mixed_radius_gregory_cube_result();
    let brep = result.brep;
    let orig_gregory_faces =
        sample_gregory_method_corner_patch_faces(&brep, &result.corner_patches);
    assert!(
        !orig_gregory_faces.is_empty(),
        "test precondition: original BRep must contain Gregory-method corner patch faces"
    );

    let doc = StepDocument {
        products: vec![StepProductNode {
            name: "GregoryFilletedCube".into(),
            id: "GFC001".into(),
            children: vec![],
            shapes: vec![Arc::new(ShapeRep {
                kind: ShapeRepKind::AdvancedBrep,
                solids: vec![brep.clone()],
                context: RepresentationContext::default(),
                non_brep_items: vec![],
            })],
            metadata: ProductMetadata::default(),
        }],
        roots: vec![0],
        units: UnitContext::default(),
        opaque_bag: OpaqueBag::default(),
    };
    let bytes = export_step(&doc, &StepExportOptions::default()).unwrap();
    // Use healing=None so the re-imported NURBS is byte-for-byte
    // the exported form; `Basic` healing can simplify / re-fit
    // geometry which would introduce additional deviation beyond
    // the rational→NURBS export path that item #2 targets.
    let import_opts = popkern_io::StepImportOptions {
        healing: popkern_io::step::healing::HealingLevel::None,
        ..Default::default()
    };
    let reimport = import_step(&bytes, &import_opts).unwrap();
    assert!(
        reimport.diagnostics.errors().next().is_none(),
        "STEP reimport must not produce errors: {:?}",
        reimport.diagnostics.entries()
    );
    let reimport_brep = &reimport.document.products[0].shapes[0].solids[0];
    let imported_nurbs: Vec<_> = reimport_brep
        .faces
        .values()
        .filter_map(|face| face.surface.as_ref())
        .filter(|sof| matches!(sof.geometry, SurfaceGeometry::Nurbs(_)))
        .collect();
    assert!(
        imported_nurbs.len() >= orig_gregory_faces.len(),
        "reimported BRep must contain at least {} NURBS faces converted from Gregory faces; got {}",
        orig_gregory_faces.len(),
        imported_nurbs.len()
    );

    // STEP exports Gregory-method cap faces as NURBS with the same
    // full-range parameterisation. Match each original cap face to the
    // imported NURBS face that minimises same-UV sample error.
    let mut max_err: f64 = 0.0;
    for original_face in &orig_gregory_faces {
        let mut best_face_err = f64::INFINITY;
        for imported in &imported_nurbs {
            let mut face_err: f64 = 0.0;
            for (orig_pt, u, v) in &original_face.samples {
                let p = imported.geometry.point_at(*u, *v);
                face_err = face_err.max(p.distance_to(orig_pt));
            }
            best_face_err = best_face_err.min(face_err);
        }
        max_err = max_err.max(best_face_err);
    }

    // Precision floor breakdown (measured): `to_nurbs_exact` at full
    // range agrees with Gregory to ~1e-10, the STEP serializer uses
    // full f64 precision per CP, and import preserves rational weights.
    assert!(
        max_err < 1e-6,
        "§5.3 item #2: Gregory-method STEP round-trip same-UV position deviation {max_err} exceeds 1e-6 tolerance (corner patch faces checked: {})",
        orig_gregory_faces.len()
    );
}
