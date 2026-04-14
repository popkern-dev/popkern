//! Tests for `g2_sweep_from_contact_solutions` — the skin-into-NurbsSurface
//! companion to the cross-section builder. Split out of
//! `g2_cross_section_tests.rs` to keep both files under the 500-line cap.

#[cfg(test)]
use crate::contact::{
    BallSide, ContactOptions, ContactRequest, ContactSolution, solve_rolling_ball_contact,
};
use crate::g2_cross_section::g2_sweep_from_contact_solutions;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};
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
