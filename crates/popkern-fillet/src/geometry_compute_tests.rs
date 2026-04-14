use super::*;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn make_cube_and_analysis() -> (BRep, analysis::EdgeAnalysis) {
    let brep = BRep::unit_cube(tol()).unwrap();
    let edge = brep.edges.keys().next().unwrap();
    let mut analyses = analysis::analyze_edges(&brep, &[edge], false).unwrap();
    (brep, analyses.remove(0))
}

#[test]
fn surface_param_range_plane_is_infinite_clamped() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let ((u_lo, u_hi), (v_lo, v_hi)) = surface_param_range(&plane);
    // Plane has infinite parameter range, so we clamp to ±10.
    assert_eq!(u_lo, -10.0);
    assert_eq!(u_hi, 10.0);
    assert_eq!(v_lo, -10.0);
    assert_eq!(v_hi, 10.0);
}

#[test]
fn surface_param_range_cylinder_is_bounded() {
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    let ((u_lo, u_hi), (v_lo, v_hi)) = surface_param_range(&cyl);
    // Cylinder: u (angle) in [0, 2π], v is infinite → clamped.
    assert!(u_lo.is_finite());
    assert!(u_hi.is_finite());
    // v axis is infinite for a cylinder, so we clamp to ±10.
    assert_eq!(v_lo, -10.0);
    assert_eq!(v_hi, 10.0);
}

#[test]
fn rebuild_analysis_preserves_handles() {
    let (brep, analysis) = make_cube_and_analysis();
    let rebuilt = rebuild_analysis(&analysis, &brep).unwrap();
    assert_eq!(rebuilt.edge, analysis.edge);
    assert_eq!(rebuilt.face_1, analysis.face_1);
    assert_eq!(rebuilt.face_2, analysis.face_2);
    assert_eq!(rebuilt.convexity, analysis.convexity);
    assert_eq!(rebuilt.vertex_start, analysis.vertex_start);
    assert_eq!(rebuilt.vertex_end, analysis.vertex_end);
}

#[test]
fn face_handle_to_index_returns_valid_index() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let faces: Vec<_> = brep.faces.keys().collect();
    for (expected_i, &face) in faces.iter().enumerate() {
        let idx = face_handle_to_index(&brep, face);
        // The index should be within the range of faces.
        assert!(idx < brep.faces.len(), "index {idx} out of range");
        // The face at that index should match.
        assert_eq!(
            brep.faces.keys().nth(idx).unwrap(),
            face,
            "index mismatch for face {expected_i}"
        );
    }
}

#[test]
fn face_handle_to_index_unknown_face_returns_zero() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let bogus = popkern_core::topology::FaceHandle::default();
    let idx = face_handle_to_index(&brep, bogus);
    assert_eq!(idx, 0, "unknown face should map to index 0");
}

#[test]
fn get_edge_midpoint_with_curve() {
    let (brep, analysis) = make_cube_and_analysis();
    let midpoint = get_edge_midpoint(&brep, &analysis).unwrap();
    // The midpoint should be within the bounding box of the cube [0,1]^3.
    assert!(midpoint.x >= -0.01 && midpoint.x <= 1.01);
    assert!(midpoint.y >= -0.01 && midpoint.y <= 1.01);
    assert!(midpoint.z >= -0.01 && midpoint.z <= 1.01);
}

#[test]
fn build_overlap_inputs_empty_geometries() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let (overlap_trims, face_trim_map) = build_overlap_inputs(&[], &brep);
    assert!(overlap_trims.is_empty());
    assert!(face_trim_map.is_empty());
}

/// Phase 2 step 4: when the caller requests G2 cross-sections on a
/// non-cone surface pair (plane-plane cube edge), `compute_fillet_geometry`
/// must route through the spine walker so that `contact_solutions` is
/// populated and `build_fillet_surface` can skin a degree-5 Bezier sweep.
/// The resulting surface must be a NURBS (not a Torus — the constant-
/// radius plane-plane analytic path), proving the G2 branch fired.
#[test]
fn compute_fillet_geometry_g2_on_cube_edge_produces_nurbs_sweep() {
    let (brep, analysis) = make_cube_and_analysis();
    let radius = crate::types::FilletRadius::Constant(0.1);
    let result = compute_fillet_geometry(&analysis, &radius, &brep, &tol(), true)
        .expect("g2 cube fillet geometry should compute");
    assert!(
        matches!(
            result.surface,
            popkern_core::geometry::surface::SurfaceGeometry::Nurbs(_)
        ),
        "g2 cube fillet should produce a NURBS sweep, not a Torus or Plane"
    );
}

#[test]
fn compute_fillet_geometry_legacy_constant_on_cube_produces_torus() {
    // Companion to the G2 test above: with `g2_cross_section = false`
    // the legacy constant-radius plane-plane analytic path runs and
    // produces an analytic Torus. This pins the behaviour the G2
    // opt-in path must NOT regress.
    let (brep, analysis) = make_cube_and_analysis();
    let radius = crate::types::FilletRadius::Constant(0.1);
    let result = compute_fillet_geometry(&analysis, &radius, &brep, &tol(), false)
        .expect("legacy cube fillet geometry should compute");
    assert!(
        matches!(
            result.surface,
            popkern_core::geometry::surface::SurfaceGeometry::Torus { .. }
        ),
        "legacy cube fillet should produce a Torus"
    );
}

#[test]
fn compute_fillet_geometry_variable_radius() {
    let (brep, analysis) = make_cube_and_analysis();
    let radius = crate::types::FilletRadius::Variable {
        start: 0.08,
        end: 0.12,
    };
    let result = compute_fillet_geometry(&analysis, &radius, &brep, &tol(), false);
    assert!(
        result.is_ok(),
        "variable radius fillet geometry should compute: {:?}",
        result.err()
    );
}

#[test]
fn compute_fillet_geometry_by_points_two_points_exercises_path() {
    // Two-point ByPoints exercises the ByPoints branch in compute_fillet_geometry.
    // The result may fail (degenerate curve), but the branch is reached.
    let (brep, analysis) = make_cube_and_analysis();
    let radius = crate::types::FilletRadius::ByPoints(vec![(0.0, 0.1), (1.0, 0.1)]);
    // Either Ok or Err is acceptable — we just confirm no panic and that the
    // ByPoints code branch fires (particularly r_start/r_end extraction).
    let _ = compute_fillet_geometry(&analysis, &radius, &brep, &tol(), false);
}

#[test]
fn compute_fillet_geometry_by_points_single_point_exercises_path() {
    let (brep, analysis) = make_cube_and_analysis();
    let radius = crate::types::FilletRadius::ByPoints(vec![(0.5, 0.1)]);
    // Single-point ByPoints uses pts[0].1 for both r_start and r_end.
    let _ = compute_fillet_geometry(&analysis, &radius, &brep, &tol(), false);
}

#[test]
fn compute_fillet_geometry_by_points_empty_fails() {
    let (brep, analysis) = make_cube_and_analysis();
    let radius = crate::types::FilletRadius::ByPoints(vec![]);
    let result = compute_fillet_geometry(&analysis, &radius, &brep, &tol(), false);
    assert!(result.is_err(), "empty ByPoints should fail");
}

#[test]
fn compute_chamfer_geometry_equal() {
    let (brep, analysis) = make_cube_and_analysis();
    let kind = crate::types::ChamferKind::Equal(0.1);
    let result = compute_chamfer_geometry(&analysis, &kind, &brep, &tol());
    assert!(
        result.is_ok(),
        "equal chamfer geometry should compute: {:?}",
        result.err()
    );
}

#[test]
fn compute_chamfer_geometry_two_distance() {
    let (brep, analysis) = make_cube_and_analysis();
    let kind = crate::types::ChamferKind::TwoDistance {
        face_1: analysis.face_1,
        dist_1: 0.1,
        dist_2: 0.2,
    };
    let result = compute_chamfer_geometry(&analysis, &kind, &brep, &tol());
    assert!(
        result.is_ok(),
        "two-distance chamfer geometry should compute: {:?}",
        result.err()
    );
}

#[test]
fn compute_chamfer_geometry_distance_angle() {
    let (brep, analysis) = make_cube_and_analysis();
    let kind = crate::types::ChamferKind::DistanceAngle {
        face: analysis.face_1,
        distance: 0.1,
        angle: std::f64::consts::FRAC_PI_4,
    };
    let result = compute_chamfer_geometry(&analysis, &kind, &brep, &tol());
    assert!(
        result.is_ok(),
        "distance-angle chamfer geometry should compute: {:?}",
        result.err()
    );
}
