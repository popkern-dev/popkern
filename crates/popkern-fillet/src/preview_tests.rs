use super::*;
use crate::types::{ChamferKind, FilletRadius};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn cube_and_edge() -> (BRep, popkern_core::topology::EdgeHandle) {
    let brep = BRep::unit_cube(tol()).unwrap();
    let edge = brep.edges.keys().next().unwrap();
    (brep, edge)
}

#[test]
fn preview_fillet_single_edge_returns_one_entry() {
    let (brep, edge) = cube_and_edge();
    let pre_face_count = brep.faces.len();

    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let preview = preview_fillet(&brep, &specs, tol()).unwrap();

    assert_eq!(preview.per_edge.len(), 1);
    assert!(preview.all_succeeded());
    // Input BRep is unchanged — preview is non-destructive.
    assert_eq!(brep.faces.len(), pre_face_count);
}

#[test]
fn preview_fillet_returns_blend_surface_and_trim_curves() {
    let (brep, edge) = cube_and_edge();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let preview = preview_fillet(&brep, &specs, tol()).unwrap();
    let pe = &preview.per_edge[0];
    assert!((pe.effective_radius - 0.1).abs() < 1e-12);
    assert_eq!(pe.edge, edge);
    // Trim curves should have non-zero extent.
    use popkern_core::geometry::curve::CurveEval;
    let (t0_1, t1_1) = pe.trim_curve_1.parameter_range();
    assert!(t1_1 > t0_1);
}

#[test]
fn preview_fillet_rejects_nan_radius() {
    let (brep, edge) = cube_and_edge();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(f64::NAN),
    }];
    assert!(matches!(
        preview_fillet(&brep, &specs, tol()),
        Err(BlendError::InvalidParameter { .. })
    ));
}

#[test]
fn preview_fillet_empty_specs_rejected() {
    let brep = BRep::unit_cube(tol()).unwrap();
    assert!(matches!(
        preview_fillet(&brep, &[], tol()),
        Err(BlendError::InvalidParameter { .. })
    ));
}

#[test]
fn preview_chamfer_returns_ruled_surface_preview() {
    let (brep, edge) = cube_and_edge();
    let specs = vec![ChamferSpec {
        edge,
        kind: ChamferKind::Equal(0.1),
    }];
    let preview = preview_chamfer(&brep, &specs, tol()).unwrap();
    assert_eq!(preview.per_edge.len(), 1);
    assert!(preview.all_succeeded());
    assert!((preview.per_edge[0].effective_radius - 0.1).abs() < 1e-12);
}

#[test]
fn preview_chamfer_rejects_negative_distance() {
    let (brep, edge) = cube_and_edge();
    let specs = vec![ChamferSpec {
        edge,
        kind: ChamferKind::Equal(-0.1),
    }];
    assert!(matches!(
        preview_chamfer(&brep, &specs, tol()),
        Err(BlendError::InvalidParameter { .. })
    ));
}

#[test]
fn preview_is_pure_and_does_not_mutate_input_brep() {
    let (brep, edge) = cube_and_edge();
    let snapshot_faces = brep.faces.len();
    let snapshot_edges = brep.edges.len();
    let snapshot_vertices = brep.vertices.len();

    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let _ = preview_fillet(&brep, &specs, tol()).unwrap();

    assert_eq!(brep.faces.len(), snapshot_faces);
    assert_eq!(brep.edges.len(), snapshot_edges);
    assert_eq!(brep.vertices.len(), snapshot_vertices);
}

/// §5.3 C7 — full-cube preview reports one corner patch per
/// vertex with ≥2 filleted edges (8 on a unit cube) and populates
/// classification + method + continuity without mutating the BRep.
#[test]
fn preview_corner_patches_on_full_cube_reports_one_per_corner() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let pre_face_count = cube.faces.len();
    let specs: Vec<FilletSpec> = cube
        .edges
        .keys()
        .map(|e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.1),
        })
        .collect();

    let previews = preview_corner_patches(&cube, &specs, tol()).unwrap();
    assert_eq!(previews.len(), 8, "cube has 8 multi-fillet vertices");
    for p in &previews {
        assert_eq!(p.n_sides, 3, "cube corner is trihedral");
        assert_eq!(
            p.classification,
            crate::vertex_classifier::VertexBlendClassification::SameSideConvex,
        );
    }
    // Non-destructive.
    assert_eq!(cube.faces.len(), pre_face_count);
}

/// §5.3 C7 — a single-edge fillet spec produces no corner
/// patches (no vertex has ≥2 filleted edges).
#[test]
fn preview_corner_patches_empty_for_single_edge_spec() {
    let (brep, edge) = cube_and_edge();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let previews = preview_corner_patches(&brep, &specs, tol()).unwrap();
    assert!(previews.is_empty());
}

#[test]
fn preview_all_succeeded_false_when_spec_fails() {
    // Craft a spec that would fail at the per-edge geometry phase:
    // an obviously-too-large radius on a curved surface that DOES
    // have a curvature bound. We can't easily do this with unit
    // cube edges (plane-plane has no bound), so this test just
    // verifies the `all_succeeded` helper behavior with a constructed
    // PreviewGeometry that has a failure.
    let pg = PreviewGeometry {
        per_edge: vec![],
        failed: vec![(
            0,
            crate::error::BlendFailure::InvalidTopology {
                reason: "test".to_string(),
                location: None,
            },
        )],
    };
    assert!(!pg.all_succeeded());
}
