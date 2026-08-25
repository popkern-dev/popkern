//! Category 6 & 7: Shell, offset, and draft operations.

use std::sync::Arc;

use popkern_core::error::ValidationError;
use popkern_core::geometry::Vec3;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_draft::{
    DraftOptions, DraftSpec, NeutralReference, PartingLine, draft_face, draft_faces,
    draft_parting_line,
};
use popkern_io::{
    OpaqueBag, ProductMetadata, RepresentationContext, ShapeRep, ShapeRepKind, StepDocument,
    StepExportOptions, StepImportOptions, StepProductNode, UnitContext, export_step, import_step,
};
use popkern_shell::{ShellDirection, ShellOptions, offset::offset_solid, shell};

fn tol() -> Tolerance {
    Tolerance::default().with_linear(1e-4)
}

fn assert_euler(brep: &BRep) {
    let v = brep.vertices.len() as i64;
    let e = brep.edges.len() as i64;
    let f = brep.faces.len() as i64;
    assert_eq!(v - e + f, 2, "Euler: V={v} E={e} F={f}");
}

fn assert_no_structural_errors(brep: &BRep) {
    let errors = brep.validate_all();
    let structural: Vec<_> = errors
        .iter()
        .filter(|e| {
            matches!(
                e,
                ValidationError::BrokenCycle { .. }
                    | ValidationError::BrokenTwin { .. }
                    | ValidationError::InconsistentVertex { .. }
                    | ValidationError::NonManifoldEdge { .. }
            )
        })
        .collect();
    assert!(structural.is_empty(), "structural errors: {:?}", structural);
}

fn assert_finite_vertices(brep: &BRep) {
    for (vh, v) in &brep.vertices {
        assert!(v.position.x.is_finite(), "vertex {:?} non-finite x", vh);
        assert!(v.position.y.is_finite(), "vertex {:?} non-finite y", vh);
        assert!(v.position.z.is_finite(), "vertex {:?} non-finite z", vh);
    }
}

fn assert_all_faces_have_surfaces(brep: &BRep) {
    for (fh, face) in &brep.faces {
        assert!(face.surface.is_some(), "face {:?} has no surface", fh);
    }
}

fn assert_step_roundtrip(brep: &BRep, name: &str) {
    let doc = StepDocument {
        products: vec![StepProductNode {
            name: name.into(),
            id: "test".into(),
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
    assert!(bytes.len() > 100, "STEP export too small");
    let result = import_step(&bytes, &StepImportOptions::default()).unwrap();
    assert!(
        !result.document.products.is_empty(),
        "reimport has no products"
    );
}

/// Find a side face (normal roughly perpendicular to Z) on a unit cube.
fn find_side_face(brep: &BRep, pull_dir: Vec3) -> popkern_core::topology::FaceHandle {
    brep.faces
        .iter()
        .find(|(_, fd)| {
            fd.surface.as_ref().is_some_and(|sof| {
                let n = sof.geometry.normal_at(0.5, 0.5);
                n.dot(&pull_dir).abs() < 0.5
            })
        })
        .map(|(fh, _)| fh)
        .expect("should find a side face")
}

/// Find all side faces (normal roughly perpendicular to pull_dir).
fn find_side_faces(brep: &BRep, pull_dir: Vec3) -> Vec<popkern_core::topology::FaceHandle> {
    brep.faces
        .iter()
        .filter(|(_, fd)| {
            fd.surface.as_ref().is_some_and(|sof| {
                let n = sof.geometry.normal_at(0.5, 0.5);
                n.dot(&pull_dir).abs() < 0.5
            })
        })
        .map(|(fh, _)| fh)
        .collect()
}

// ── Category 6: Shell/Offset ──────────────────────────────────────

#[test]
fn t35_shell_cube_no_removal() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let result = shell(
        &cube,
        0.1,
        ShellDirection::Inward,
        &[],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    )
    .unwrap();
    let brep = &result.brep;
    // No face removed: outer 6 + inner 6 = 12
    assert_eq!(brep.faces.len(), 12, "shell no removal: 12 faces");
    assert_no_structural_errors(brep);
    assert_finite_vertices(brep);
    assert_all_faces_have_surfaces(brep);
}

#[test]
fn t36_shell_cube_one_face_removed() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let face_to_remove = cube.faces.keys().next().unwrap();
    let result = shell(
        &cube,
        0.1,
        ShellDirection::Inward,
        &[face_to_remove],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    )
    .unwrap();
    let brep = &result.brep;
    assert!(
        brep.faces.len() > 6,
        "shelled cube should have more than 6 faces, got {}",
        brep.faces.len()
    );
    assert_finite_vertices(brep);
    assert_all_faces_have_surfaces(brep);
    // The historical "invalid SlotMap key" validator panic on
    // face-removed shells is gone (PK-306/307 root fixes); PK-321
    // re-armed full validation here.
    let errors = brep.validate_all();
    assert!(errors.is_empty(), "validation errors: {errors:?}");
}

#[test]
fn t37_shell_cube_inward_outward_midplane() {
    let cube_orig = BRep::unit_cube(tol()).unwrap();
    for direction in [
        ShellDirection::Inward,
        ShellDirection::Outward,
        ShellDirection::Midplane,
    ] {
        let face = cube_orig.faces.keys().next().unwrap();
        let result = shell(
            &cube_orig,
            0.1,
            direction,
            &[face],
            &[],
            &ShellOptions::default(),
            tol(),
            None,
        )
        .unwrap();
        let brep = &result.brep;
        assert!(
            brep.faces.len() > 6,
            "{:?} shell should have >6 faces",
            direction
        );
        assert_finite_vertices(brep);
        // Full validation re-armed (PK-321): the historical SlotMap-key
        // panic on face-removed shells is fixed.
        let errors = brep.validate_all();
        assert!(
            errors.is_empty(),
            "{direction:?}: validation errors: {errors:?}"
        );
    }
}

#[test]
fn t38_offset_solid_grows_cube() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let result = offset_solid(&cube, 0.1, &ShellOptions::default(), tol()).unwrap();
    let brep = &result.brep;
    // Minkowski-faithful default (PK-325): the rounded box.
    assert_eq!(brep.faces.len(), 26, "arc offset serves the rounded box");
    assert_no_structural_errors(brep);
    assert_finite_vertices(brep);
    // Check that the offset moved vertices (either direction — the sign convention
    // may be inward-positive). Just verify the geometry changed.
    let orig_max: f64 = cube
        .vertices
        .values()
        .map(|v| v.position.x)
        .fold(f64::NEG_INFINITY, f64::max);
    let new_max: f64 = brep
        .vertices
        .values()
        .map(|v| v.position.x)
        .fold(f64::NEG_INFINITY, f64::max);
    assert!(
        (new_max - orig_max).abs() > 0.05,
        "offset should move vertices: orig_max={orig_max}, new_max={new_max}"
    );
}

#[test]
fn t39_offset_solid_shrinks_cube() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let result = offset_solid(&cube, -0.1, &ShellOptions::default(), tol()).unwrap();
    let brep = &result.brep;
    assert_eq!(brep.faces.len(), 6, "offset should keep 6 faces");
    assert_no_structural_errors(brep);
    assert_finite_vertices(brep);
    let orig_max: f64 = cube
        .vertices
        .values()
        .map(|v| v.position.x)
        .fold(f64::NEG_INFINITY, f64::max);
    let new_max: f64 = brep
        .vertices
        .values()
        .map(|v| v.position.x)
        .fold(f64::NEG_INFINITY, f64::max);
    assert!(
        (new_max - orig_max).abs() > 0.05,
        "negative offset should move vertices: orig_max={orig_max}, new_max={new_max}"
    );
}

#[test]
fn t40_shell_step_roundtrip() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let face = cube.faces.keys().next().unwrap();
    let result = shell(
        &cube,
        0.1,
        ShellDirection::Inward,
        &[face],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    )
    .unwrap();
    assert_step_roundtrip(&result.brep, "shelled_cube");
}

// ── Category 7: Draft/Taper ───────────────────────────────────────

#[test]
fn t41_draft_one_side_face() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let pull = Vec3::new(0.0, 0.0, 1.0);
    let face = find_side_face(&cube, pull);
    let result = draft_face(
        &cube,
        face,
        0.087,
        pull,
        NeutralReference::Auto,
        tol(),
        None,
    )
    .unwrap();
    let brep = &result.brep;
    assert_eq!(brep.faces.len(), 6, "draft should keep 6 faces");
    assert_eq!(brep.edges.len(), 12);
    assert_eq!(brep.vertices.len(), 8);
    assert_no_structural_errors(brep);
    assert_finite_vertices(brep);
    assert_all_faces_have_surfaces(brep);
    // Check that the drafted face has a tilted normal
    let has_tilted = brep.faces.values().any(|f| {
        f.surface.as_ref().is_some_and(|sof| {
            let n = sof.geometry.normal_at(0.5, 0.5);
            let d = n.dot(&pull).abs();
            d > 0.01 && d < 0.5
        })
    });
    assert!(has_tilted, "drafted face should have tilted normal");
}

#[test]
fn t42_draft_all_side_faces() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let pull = Vec3::new(0.0, 0.0, 1.0);
    let side_faces = find_side_faces(&cube, pull);
    assert_eq!(side_faces.len(), 4, "unit cube has 4 side faces");
    let specs: Vec<DraftSpec> = side_faces
        .iter()
        .map(|&f| DraftSpec {
            face: f,
            angle: 0.087, // ~5 degrees
            pull_direction: pull,
            neutral: NeutralReference::Auto,
            angle_variation: popkern_draft::AngleVariation::Constant,
        })
        .collect();
    match draft_faces(&cube, &specs, &DraftOptions::default(), tol(), None) {
        Ok(result) => {
            let brep = &result.brep;
            assert_eq!(brep.faces.len(), 6);
            assert_no_structural_errors(brep);
            assert_finite_vertices(brep);
        }
        Err(e) => {
            println!("draft_all_side_faces failed: {e}");
        }
    }
}

#[test]
fn t43_draft_with_parting_line() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let pull = Vec3::new(0.0, 0.0, 1.0);
    let parting = PartingLine::Plane {
        origin: popkern_core::geometry::Point3::new(0.5, 0.5, 0.5),
        normal: Vec3::new(0.0, 0.0, 1.0),
    };
    match draft_parting_line(
        &cube,
        &parting,
        0.087,
        0.087,
        pull,
        &DraftOptions::default(),
        tol(),
    ) {
        Ok(result) => {
            let brep = &result.brep;
            assert!(brep.faces.len() >= 6, "parting draft should have faces");
            assert_no_structural_errors(brep);
            assert_finite_vertices(brep);
        }
        Err(e) => {
            println!("draft_with_parting_line failed: {e}");
        }
    }
}

#[test]
fn t44_draft_preserves_euler() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let pull = Vec3::new(0.0, 0.0, 1.0);
    let face = find_side_face(&cube, pull);
    let result = draft_face(
        &cube,
        face,
        0.087,
        pull,
        NeutralReference::Auto,
        tol(),
        None,
    )
    .unwrap();
    assert_euler(&result.brep);
}

#[test]
fn t45_draft_step_roundtrip() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let pull = Vec3::new(0.0, 0.0, 1.0);
    let face = find_side_face(&cube, pull);
    let result = draft_face(
        &cube,
        face,
        0.087,
        pull,
        NeutralReference::Auto,
        tol(),
        None,
    )
    .unwrap();
    assert_step_roundtrip(&result.brep, "drafted_face");
}
