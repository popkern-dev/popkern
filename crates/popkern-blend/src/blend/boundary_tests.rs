//! Tests for `blend::boundary`. Extracted so the implementation
//! file stays under the 500-line structural cap.

use super::boundary::{BoundarySegment, compute_blend_boundary};
use crate::fillet_surface::FilletSurfaceResult;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn test_tolerance() -> Tolerance {
    Tolerance::default()
        .with_linear(1e-3)
        .with_angular(1e-3)
        .with_parametric(1e-6)
}

fn make_fillet_result(
    start1: Point3,
    end1: Point3,
    start2: Point3,
    end2: Point3,
) -> FilletSurfaceResult {
    let trim_1 =
        crate::center_curve::fit_nurbs_through_points(&[start1, end1]).expect("fit should succeed");
    let trim_2 =
        crate::center_curve::fit_nurbs_through_points(&[start2, end2]).expect("fit should succeed");
    FilletSurfaceResult {
        surface: SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        },
        trim_curve_1: trim_1,
        trim_curve_2: trim_2,
        continuity_downgrade: None,
    }
}

#[test]
fn blend_boundary_two_fillets() {
    let tol = test_tolerance();
    let brep = BRep::unit_cube(tol).expect("cube should build");
    let vertex = brep.vertices.keys().next().expect("cube has vertices");
    let vertex_pos = brep.vertices.get(vertex).expect("vertex exists").position;

    let f1 = make_fillet_result(
        vertex_pos,
        Point3::new(vertex_pos.x + 1.0, vertex_pos.y, vertex_pos.z),
        Point3::new(vertex_pos.x, vertex_pos.y + 0.1, vertex_pos.z),
        Point3::new(vertex_pos.x + 1.0, vertex_pos.y + 0.1, vertex_pos.z),
    );
    let f2 = make_fillet_result(
        vertex_pos,
        Point3::new(vertex_pos.x, vertex_pos.y + 1.0, vertex_pos.z),
        Point3::new(vertex_pos.x + 0.1, vertex_pos.y, vertex_pos.z),
        Point3::new(vertex_pos.x + 0.1, vertex_pos.y + 1.0, vertex_pos.z),
    );

    let fillets = vec![f1, f2];
    let result = compute_blend_boundary(vertex, &[0, 1], &fillets, &brep, &tol);
    assert!(result.is_ok(), "blend boundary failed: {:?}", result.err());

    let boundary = result.expect("already checked");
    // §5.3: `compute_blend_boundary` emits alternating
    // FilletEnd / FaceEdge segments on planar-face vertices. The
    // test's synthetic plane fillets may or may not trigger the
    // bridge path depending on the cube-face projection, so the
    // assertion is on FilletEnd count (stable) rather than total
    // segment count.
    let n_fe = boundary
        .segments
        .iter()
        .filter(|s| matches!(s, BoundarySegment::FilletEnd { .. }))
        .count();
    assert_eq!(
        n_fe, 2,
        "expected 2 FilletEnd segments regardless of face-edge bridges",
    );
    assert_eq!(boundary.vertex, vertex);
}

#[test]
fn blend_boundary_empty_fillets() {
    let tol = test_tolerance();
    let brep = BRep::unit_cube(tol).expect("cube should build");
    let vertex = brep.vertices.keys().next().expect("cube has vertices");
    let result = compute_blend_boundary(vertex, &[], &[], &brep, &tol);
    assert!(result.is_err(), "should fail with no fillet indices");
}

#[test]
fn blend_boundary_invalid_fillet_index() {
    let tol = test_tolerance();
    let brep = BRep::unit_cube(tol).expect("cube should build");
    let vertex = brep.vertices.keys().next().expect("cube has vertices");
    let result = compute_blend_boundary(vertex, &[99], &[], &brep, &tol);
    assert!(result.is_err(), "should fail with out-of-range index");
}

/// PK-548/PK-549: the closest-to-vertex station (the ring push-back's
/// and the default fillet-end station) is the exact foot of the
/// perpendicular. It is independent of the rail's direction, its
/// extension past the vertex, and the sample count. The quantized
/// search it replaces put a reversed rail's station 0.0066 off its
/// neighbours'.
#[test]
fn the_closest_station_is_exact_and_direction_independent() {
    use popkern_core::geometry::NurbsCurve;
    let line = |a: Point3, b: Point3| NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![a, b],
        weights: vec![1.0, 1.0],
    };
    let vertex = Point3::new(0.0, 0.0, 0.8);
    // A rail parallel to the edge, offset, overshooting the vertex
    // differently in each direction.
    let (a, b) = (
        Point3::new(-0.848, 0.03, 1.4573),
        Point3::new(1.848, 0.03, -0.6997),
    );
    let direction = b - a;
    let foot_t = (vertex - a).dot(&direction) / direction.length_squared();
    let foot = a + direction * foot_t;
    for samples in [16usize, 64, 128, 512] {
        for rail in [line(a, b), line(b, a)] {
            let t = super::boundary_hex::closest_param_on_nurbs_curve(&rail, &vertex, samples);
            let p = popkern_core::geometry::curve::CurveEval::point_at(&rail, t);
            assert!(
                p.distance_to(&foot) < 1e-12,
                "samples {samples}: {p:?} vs {foot:?}"
            );
        }
    }
}
