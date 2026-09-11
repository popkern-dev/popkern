//! Overflow of `split_extra_tests.rs` (500-line structural budget).

use super::*;
use popkern_core::{
    geometry::{Point3, nurbs_curve::NurbsCurve},
    tolerance::Tolerance,
    topology::BRep,
};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Face with pcurves but where the UV curve doesn't cross the boundary
/// should be skipped gracefully (0 crossings).
#[test]
fn uv_curve_outside_boundary_skips_face() {
    use popkern_core::geometry::curve_on_surface::CurveOnSurface;

    let cube = BRep::unit_cube(tol()).expect("cube");
    let face_h = cube.faces.keys().next().expect("has faces");

    // Curve entirely outside the cube's UV domain: u,v in [10,11].
    // Use 3D coords far off-plane so the planar fallback also rejects it.
    let far_uv = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(10.0, 10.0, 0.0), Point3::new(11.0, 11.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let far_3d = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(10.0, 10.0, 5.0), Point3::new(11.0, 11.0, 5.0)],
        weights: vec![1.0, 1.0],
    };
    let curve = CurveOnSurface {
        curve_3d: far_3d,
        curve_uv1: far_uv.clone(),
        curve_uv2: far_uv,
        is_closed: false,
    };
    let record = IntersectionRecord {
        residue_owner_a: None,
        face_a: face_h,
        face_b: popkern_core::topology::FaceHandle::default(),
        curve,
        is_tangent: false,
        is_coincident: false,
        start_boundary: Default::default(),
        end_boundary: Default::default(),
    };

    let brep_b = BRep::new(tol());
    let faces_before = cube.faces.len();
    let result = split_along_curves(cube, brep_b, &[record], &tol())
        .expect("should succeed, skipping the face");

    assert_eq!(
        result.a.faces.len(),
        faces_before,
        "face should not be split when UV curve is outside boundary"
    );
}
