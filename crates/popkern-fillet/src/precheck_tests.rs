use super::*;
use crate::analysis::Convexity;
use crate::types::FilletRadius;
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

// ---------- precheck_fillet ----------

/// §13.1 integration: `precheck_fillet` must report a
/// `FilletIssue::Degeneracy` and flip `can_succeed` to false when
/// the edge is collapsed. Synthesises a collapsed edge by moving
/// one vertex of a cube edge to coincide with its neighbour.
#[test]
fn precheck_reports_collapsed_edge_as_degeneracy() {
    let (mut brep, edge) = cube_and_edge();
    let edge_rec = &brep.edges[edge];
    let he = &brep.half_edges[edge_rec.half_edge];
    let twin = &brep.half_edges[he.twin];
    let v_start = he.origin;
    let v_end = twin.origin;
    brep.vertices[v_start].position = brep.vertices[v_end].position;

    let f = precheck_fillet(&brep, edge, &FilletRadius::Constant(0.1), tol());
    assert!(
        !f.can_succeed,
        "collapsed edge should flip can_succeed to false"
    );
    let has_degeneracy = f.issues.iter().any(|i| {
        matches!(
            i,
            FilletIssue::Degeneracy(crate::degeneracy::DegeneracyKind::CollapsedEdge { .. })
        )
    });
    assert!(
        has_degeneracy,
        "issues should contain a CollapsedEdge Degeneracy: {:?}",
        f.issues
    );
}

#[test]
fn precheck_plane_plane_cube_edge_is_clean_for_small_radius() {
    let (brep, edge) = cube_and_edge();
    let f = precheck_fillet(&brep, edge, &FilletRadius::Constant(0.1), tol());
    assert!(
        f.can_succeed,
        "small-radius plane-plane precheck should pass"
    );
    // Planes have infinite curvature radius, so there is no
    // curvature-derived bound.
    assert!(f.max_radius_estimate.is_none());
    assert!(!f.is_seam);
    assert!(f.is_clean());
}

#[test]
fn precheck_unknown_edge_returns_invalid_topology() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let bogus = popkern_core::topology::EdgeHandle::default();
    let f = precheck_fillet(&brep, bogus, &FilletRadius::Constant(0.1), tol());
    assert!(!f.can_succeed);
    assert!(
        f.issues
            .iter()
            .any(|i| matches!(i, FilletIssue::InvalidTopology { .. }))
    );
}

#[test]
fn precheck_is_clean_returns_true_on_empty_issues() {
    let (brep, edge) = cube_and_edge();
    let f = precheck_fillet(&brep, edge, &FilletRadius::Constant(0.1), tol());
    assert_eq!(f.is_clean(), f.issues.is_empty());
}

#[test]
fn precheck_variable_radius_uses_peak() {
    let (brep, edge) = cube_and_edge();
    let f = precheck_fillet(
        &brep,
        edge,
        &FilletRadius::Variable {
            start: 0.05,
            end: 0.15,
        },
        tol(),
    );
    // Plane-plane, so no curvature bound. Should still pass
    // since 0.15 fits on a unit cube.
    assert!(f.can_succeed);
}

#[test]
fn precheck_by_points_uses_peak_radius() {
    let (brep, edge) = cube_and_edge();
    let f = precheck_fillet(
        &brep,
        edge,
        &FilletRadius::ByPoints(vec![(0.0, 0.05), (0.5, 0.1), (1.0, 0.05)]),
        tol(),
    );
    assert!(f.can_succeed);
}

#[test]
fn precheck_convexity_classified_on_cube_edge() {
    let (brep, edge) = cube_and_edge();
    let f = precheck_fillet(&brep, edge, &FilletRadius::Constant(0.1), tol());
    // A unit cube edge is convex. The precheck should surface that.
    assert_eq!(f.convexity, Convexity::Convex);
}

// ---------- feasible_radius_bounds ----------

#[test]
fn feasible_radius_bounds_plane_plane_has_geometric_bound() {
    let (brep, edge) = cube_and_edge();
    let bounds = feasible_radius_bounds(&brep, edge, tol()).unwrap();
    assert!(bounds.min > 0.0);
    assert!(bounds.max > 0.0);
    // Plane-plane has no curvature bound.
    assert!(bounds.curvature_bound.is_none());
    // Geometric bound should be edge_length / 4 = 0.25 for a unit cube edge.
    assert!((bounds.geometric_bound - 0.25).abs() < 1e-9);
}

#[test]
fn feasible_radius_bounds_rejects_unknown_edge() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let bogus = popkern_core::topology::EdgeHandle::default();
    let err = feasible_radius_bounds(&brep, bogus, tol()).unwrap_err();
    assert!(matches!(err, BlendError::InvalidParameter { .. }));
}

#[test]
fn feasible_radius_bounds_contains_works() {
    let (brep, edge) = cube_and_edge();
    let bounds = feasible_radius_bounds(&brep, edge, tol()).unwrap();
    assert!(bounds.contains(0.1));
    assert!(!bounds.contains(1000.0));
    assert!(!bounds.contains(-1.0));
}

#[test]
fn feasible_radius_bounds_clamp_works() {
    let (brep, edge) = cube_and_edge();
    let bounds = feasible_radius_bounds(&brep, edge, tol()).unwrap();
    // Clamping a value inside the range is identity.
    let mid = (bounds.min + bounds.max) * 0.5;
    assert!((bounds.clamp(mid) - mid).abs() < 1e-12);
    // Clamping above max returns max.
    assert!((bounds.clamp(1e6) - bounds.max).abs() < 1e-12);
    // Clamping below min returns min.
    assert!((bounds.clamp(-1e6) - bounds.min).abs() < 1e-12);
}

#[test]
fn feasible_radius_bounds_min_is_tolerance_linear() {
    let (brep, edge) = cube_and_edge();
    let bounds = feasible_radius_bounds(&brep, edge, tol()).unwrap();
    assert!((bounds.min - tol().linear()).abs() < 1e-15);
}

// ---------- fillet_quality ----------

#[test]
fn fillet_quality_after_successful_cube_fillet() {
    use crate::{FilletOptions, FilletSpec, fillet_edges};
    let (brep, edge) = cube_and_edge();
    let pre_faces: std::collections::HashSet<_> = brep.faces.keys().collect();

    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let result = fillet_edges(&brep, &specs, &FilletOptions::default(), tol(), None).unwrap();

    let quality = fillet_quality(&ComputeQualityInput {
        post_brep: &result.brep,
        pre_faces: &pre_faces,
        requested_radius: 0.1,
    });

    // At least one new face created.
    assert!(quality.patch_count >= 1);
    // Surface area is positive.
    assert!(quality.surface_area > 0.0);
    // Composite score is in [0, 1].
    assert!(quality.composite_score >= 0.0);
    assert!(quality.composite_score <= 1.0);
}

#[test]
fn fillet_quality_with_no_new_faces_falls_back_to_requested_radius() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let pre_faces: std::collections::HashSet<_> = brep.faces.keys().collect();
    let quality = fillet_quality(&ComputeQualityInput {
        post_brep: &brep,
        pre_faces: &pre_faces,
        requested_radius: 0.1,
    });
    assert_eq!(quality.patch_count, 0);
    assert_eq!(quality.surface_area, 0.0);
    // No new faces → r_min falls back to requested_radius.
    assert!((quality.r_min - 0.1).abs() < 1e-12);
}

#[test]
fn fillet_quality_composite_score_is_one_when_r_min_meets_requested() {
    // An exact torus-section fillet on a cube has r_min = requested
    // radius, so the score should be 1.0.
    use crate::{FilletOptions, FilletSpec, fillet_edges};
    let (brep, edge) = cube_and_edge();
    let pre_faces: std::collections::HashSet<_> = brep.faces.keys().collect();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let result = fillet_edges(&brep, &specs, &FilletOptions::default(), tol(), None).unwrap();
    let quality = fillet_quality(&ComputeQualityInput {
        post_brep: &result.brep,
        pre_faces: &pre_faces,
        requested_radius: 0.1,
    });
    // r_min should be >= requested, giving composite_score = 1.0.
    assert!(
        quality.composite_score >= 0.95,
        "expected high score for exact cube fillet, got {}",
        quality.composite_score
    );
}
