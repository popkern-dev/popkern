//! Tests for `g2_sweep_from_contact_solutions` — the skin-into-NurbsSurface
//! companion to the cross-section builder. Split out of
//! `g2_cross_section_tests.rs` to keep both files under the 500-line cap.

#[cfg(test)]
use crate::contact::{
    BallSide, ContactOptions, ContactRequest, ContactSolution, solve_rolling_ball_contact,
};
use crate::g2_cross_section::g2_sweep_from_contact_solutions;
use crate::walker::{WalkOptions, WalkRequest, walk_spine};
use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{NurbsSurface, Point3, Vec3};
use popkern_core::tolerance::Tolerance;

fn unit_tol() -> Tolerance {
    Tolerance::default()
}

#[test]
fn g2_sweep_skins_walker_solutions_into_nurbs_surface() {
    // Build a small set of plane-plane contact solutions (3 spine
    // samples along the y-axis), then skin them into a G2 NurbsSurface.
    // The surface should be degree-1 in u (skin direction), degree-5
    // in v (cross-section direction), with 3·6 = 18 control points.
    let plane_z0 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
    };
    let plane_x0 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        u_axis: Vec3::new(0.0, 1.0, 0.0),
        v_axis: Vec3::new(0.0, 0.0, 1.0),
        normal: Vec3::new(1.0, 0.0, 0.0),
    };
    let radius = 0.5;
    let mut solutions = Vec::new();
    for i in 0..3 {
        let y = i as f64;
        let request = ContactRequest {
            surface_0: &plane_z0,
            surface_1: &plane_x0,
            side_0: BallSide::Outside,
            side_1: BallSide::Outside,
            radius,
            initial_uv0: (radius, y),
            uv0_range: ((-10.0, 10.0), (-10.0, 10.0)),
            initial_uv1: (y, radius),
            uv1_range: ((-10.0, 10.0), (-10.0, 10.0)),
            spine_point: Point3::new(0.0, y, 0.0),
            spine_plane_normal: Vec3::new(0.0, 1.0, 0.0),
        };
        let solution =
            solve_rolling_ball_contact(&request, &ContactOptions::default(), &unit_tol(), None)
                .expect("plane-plane should converge");
        solutions.push(solution);
    }

    let surface = g2_sweep_from_contact_solutions(
        &solutions,
        &plane_z0,
        &plane_x0,
        BallSide::Outside,
        BallSide::Outside,
    )
    .expect("G2 sweep should build");

    assert_eq!(surface.degree_u, 1);
    assert_eq!(surface.degree_v, 5);
    assert_eq!(surface.num_u, 3);
    assert_eq!(surface.num_v, 6);
    assert_eq!(surface.control_points.len(), 18);
    assert_eq!(surface.weights.len(), 18);

    // The first row of control points (i=0, j=0..5) should be the
    // first cross-section's control points. Verify P_0 (the first
    // cp) equals the first solution's contact0.
    assert!(
        (surface.control_points[0] - solutions[0].contact0).length() < 1e-10,
        "first control point {:?} ≠ first sample contact0 {:?}",
        surface.control_points[0],
        solutions[0].contact0
    );
    // And the last cp of the first row should equal contact1.
    assert!(
        (surface.control_points[5] - solutions[0].contact1).length() < 1e-10,
        "first row last cp {:?} ≠ first sample contact1 {:?}",
        surface.control_points[5],
        solutions[0].contact1
    );

    // The first row (u=0) sampled at v=0 should give the first
    // solution's contact0.
    let p_at_0_0 = surface.point_at(0.0, 0.0);
    assert!(
        (p_at_0_0 - solutions[0].contact0).length() < 1e-10,
        "S(0,0) {:?} ≠ first sample contact0 {:?}",
        p_at_0_0,
        solutions[0].contact0
    );
}

/// Extract the `i`-th u-row of control points from a G2-skinned
/// NurbsSurface as a standalone degree-5 Bezier control polygon.
/// The skinned surface is linear in u and degree-5 in v, so each
/// u-row (six consecutive CPs) IS a single cross-section Bezier.
fn cross_section_row(surface: &NurbsSurface, i: usize) -> [Point3; 6] {
    let base = i * surface.num_v;
    [
        surface.control_points[base],
        surface.control_points[base + 1],
        surface.control_points[base + 2],
        surface.control_points[base + 3],
        surface.control_points[base + 4],
        surface.control_points[base + 5],
    ]
}

/// Curvature of a degree-5 Bezier at its start (t = 0), computed
/// directly from the control polygon.
///
/// For a degree-n Bezier:
/// - `P'(0) = n · (P_1 − P_0)`
/// - `P''(0) = n(n−1) · (P_2 − 2·P_1 + P_0)`
///
/// Then `κ = |P' × P''| / |P'|³`.
fn bezier5_curvature_at_start(cp: &[Point3; 6]) -> f64 {
    let p_prime = (cp[1] - cp[0]) * 5.0;
    // P''(0) = 20 · (P_2 − 2·P_1 + P_0) = 20 · ((P_2 − P_1) − (P_1 − P_0))
    let p_double_prime = ((cp[2] - cp[1]) - (cp[1] - cp[0])) * 20.0;
    let cross = p_prime.cross(&p_double_prime);
    let speed3 = p_prime.length().powi(3);
    if speed3 < 1e-15 {
        return 0.0;
    }
    cross.length() / speed3
}

/// Curvature of a degree-5 Bezier at its end (t = 1).
fn bezier5_curvature_at_end(cp: &[Point3; 6]) -> f64 {
    let p_prime = (cp[5] - cp[4]) * 5.0;
    // P''(1) = 20 · (P_5 − 2·P_4 + P_3) = 20 · ((P_5 − P_4) − (P_4 − P_3))
    let p_double_prime = ((cp[5] - cp[4]) - (cp[4] - cp[3])) * 20.0;
    let cross = p_prime.cross(&p_double_prime);
    let speed3 = p_prime.length().powi(3);
    if speed3 < 1e-15 {
        return 0.0;
    }
    cross.length() / speed3
}

/// §2.4 end-to-end acceptance: run the walker on a real plane-
/// cylinder configuration, skin the resulting contact solutions
/// into a NurbsSurface via `g2_sweep_from_contact_solutions`, then
/// extract each cross-section row from the skinned surface's
/// control points and verify every row meets the §2.4 1% curvature-
/// matching threshold at **both** endpoints.
///
/// This is the full walker → sweep → skinned-surface pipeline test
/// that `g2_cross_section_curve_curvature_meets_2_4_acceptance_on_cylinder`
/// leaves uncovered — that test builds a single cross-section from a
/// single solution; this test verifies every cross-section in the
/// *skinned* surface still has the right endpoint curvature after
/// the control points have been laid out into the surface's
/// `control_points` row-major buffer.
#[test]
fn g2_sweep_end_to_end_meets_2_4_acceptance_on_plane_cylinder_walker() {
    // Configuration: cylinder along z with radius 2, plane y = 0
    // (the xz plane). Intersection is the line x = 2, y = 0. Walker
    // spine runs along that line in z.
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 2.0,
    };
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 0.0, 1.0),
        normal: Vec3::new(0.0, 1.0, 0.0),
    };
    let spine = CurveGeometry::Line {
        origin: Point3::new(2.0, 0.0, 0.0),
        direction: Vec3::new(0.0, 0.0, 1.0),
    };
    let r = 0.2;
    let radius_at = |_t: f64| r;

    let walk_request = WalkRequest {
        surface_0: &cyl,
        surface_1: &plane,
        side_0: BallSide::Outside,
        side_1: BallSide::Outside,
        radius_at: &radius_at,
        spine: &spine,
        t_range: (0.0, 1.0),
        uv0_range: ((-3.2, 6.3), (-10.0, 10.0)),
        uv1_range: ((-10.0, 10.0), (-10.0, 10.0)),
        initial_seed: None,
        surface_0_periodicity: cyl.periodicity(),
        surface_1_periodicity: plane.periodicity(),
    };
    let walk = walk_spine(&walk_request, &WalkOptions::default(), &unit_tol())
        .expect("plane-cylinder walk should succeed");
    assert!(
        walk.samples.len() >= 3,
        "walker should produce ≥3 spine samples for a meaningful sweep"
    );

    // Skin the real walker solutions into a G2 NurbsSurface.
    let solutions: Vec<ContactSolution> = walk.samples.iter().map(|s| s.solution.clone()).collect();
    let surface = g2_sweep_from_contact_solutions(
        &solutions,
        &cyl,
        &plane,
        BallSide::Outside,
        BallSide::Outside,
    )
    .expect("plane-cylinder G2 sweep should build");

    assert_eq!(surface.degree_u, 1);
    assert_eq!(surface.degree_v, 5);
    assert_eq!(surface.num_u, solutions.len());
    assert_eq!(surface.num_v, 6);

    // §2.4 acceptance threshold: |κ_curve − κ_face| / max(κ_curve, κ_face, 1) < 0.01.
    // At the cylinder endpoint (t = 0, surface_0 = cylinder side) the
    // cross-section tangent is purely in the cylinder's angular
    // direction (see the companion test for the derivation), so
    // κ_face = 1/R = 0.5. At the plane endpoint (t = 1, surface_1 =
    // plane side) the plane has zero curvature everywhere.
    let expected_kappa_cyl = 1.0 / 2.0;
    let expected_kappa_plane = 0.0;
    let threshold = 0.01;

    for i in 0..surface.num_u {
        let row = cross_section_row(&surface, i);

        let kappa_cyl = bezier5_curvature_at_start(&row);
        let denom_cyl = kappa_cyl.max(expected_kappa_cyl).max(1.0);
        let rel_err_cyl = (kappa_cyl - expected_kappa_cyl).abs() / denom_cyl;
        assert!(
            rel_err_cyl < threshold,
            "row {i}: cylinder-side curvature {kappa_cyl:.6} \
             deviates from 1/R = {expected_kappa_cyl} by \
             {rel_err_cyl:.4} (threshold {threshold})"
        );

        let kappa_plane = bezier5_curvature_at_end(&row);
        // For the plane side the face curvature is zero, so the
        // acceptance reduces to an absolute threshold on |κ_curve|.
        assert!(
            (kappa_plane - expected_kappa_plane).abs() < threshold,
            "row {i}: plane-side curvature {kappa_plane:.6} should be ≈ 0 \
             (threshold {threshold})"
        );
    }
}

#[test]
fn g2_sweep_rejects_solutions_with_fewer_than_two_samples() {
    let solutions: Vec<ContactSolution> = Vec::new();
    let plane_z0 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
    };
    let plane_x0 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        u_axis: Vec3::new(0.0, 1.0, 0.0),
        v_axis: Vec3::new(0.0, 0.0, 1.0),
        normal: Vec3::new(1.0, 0.0, 0.0),
    };
    let result = g2_sweep_from_contact_solutions(
        &solutions,
        &plane_z0,
        &plane_x0,
        BallSide::Outside,
        BallSide::Outside,
    );
    assert!(result.is_err());
}
