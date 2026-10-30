//! Edge-in-face interference: spans, boundary exclusion, coverage.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use super::recover_edge_in_face_records;
use crate::intersect::{IntersectionRecord, compute_face_aabbs};

/// The box [x0,x1]×[y0,y1]×[z0,z1].
fn block(lo: (f64, f64, f64), hi: (f64, f64, f64), tol: Tolerance) -> BRep {
    let profile = [
        Point3::new(lo.0, lo.1, lo.2),
        Point3::new(hi.0, lo.1, lo.2),
        Point3::new(hi.0, hi.1, lo.2),
        Point3::new(lo.0, hi.1, lo.2),
    ];
    BRep::prism(&profile, Vec3::new(0.0, 0.0, 1.0), hi.2 - lo.2, tol).expect("block")
}

/// The planar face through `p` whose plane normal is ±`axis`.
fn face_at(brep: &BRep, p: Point3, axis: usize) -> FaceHandle {
    brep.faces
        .iter()
        .find(|(_, f)| match f.surface.as_ref().map(|s| &s.geometry) {
            Some(SurfaceGeometry::Plane { origin, normal, .. }) => {
                let n = [normal.x, normal.y, normal.z];
                n[axis].abs() > 0.999 && (p - *origin).dot(normal).abs() < 1e-9
            }
            _ => false,
        })
        .map(|(h, _)| h)
        .expect("face")
}

fn ends(r: &IntersectionRecord) -> (Point3, Point3) {
    let c = &r.curve.curve_3d;
    let (t0, t1) = c.parameter_range();
    (c.point_at(t0), c.point_at(t1))
}

fn runs_between(r: &IntersectionRecord, p: Point3, q: Point3) -> bool {
    let (s, e) = ends(r);
    (s.distance_to(&p) < 1e-9 && e.distance_to(&q) < 1e-9)
        || (s.distance_to(&q) < 1e-9 && e.distance_to(&p) < 1e-9)
}

fn recover(a: &BRep, b: &BRep, records: &mut Vec<IntersectionRecord>, tol: &Tolerance) {
    let (ba, bb) = (compute_face_aabbs(a), compute_face_aabbs(b));
    recover_edge_in_face_records(a, b, &ba, &bb, records, tol);
}

/// A rests on B, its bottom on B's top, the two overlapping in part:
/// A's bottom edge x = 2 crosses B's top face between y = 0.5 and 1.5,
/// and that span alone is imprinted, ending exactly where the edge
/// crosses B's boundary. B's top edges inside A's bottom face come back
/// the same way; the coplanar pair itself gets nothing.
#[test]
fn the_span_of_an_edge_inside_a_face_is_imprinted_between_its_boundary_crossings() {
    let tol = Tolerance::default();
    let a = block((0.0, 0.0, 0.0), (2.0, 2.0, 1.0), tol);
    let b = block((1.0, 0.5, -1.0), (3.0, 1.5, 0.0), tol);
    let mut records = Vec::new();
    recover(&a, &b, &mut records, &tol);

    let a_side = face_at(&a, Point3::new(2.0, 1.0, 0.5), 0);
    let a_bottom = face_at(&a, Point3::new(1.0, 1.0, 0.0), 2);
    let b_top = face_at(&b, Point3::new(2.0, 1.0, 0.0), 2);
    let b_left = face_at(&b, Point3::new(1.0, 1.0, -0.5), 0);
    let b_front = face_at(&b, Point3::new(2.0, 0.5, -0.5), 1);

    let span = |fa, fb, p, q| {
        records
            .iter()
            .filter(|r| r.face_a == fa && r.face_b == fb && runs_between(r, p, q))
            .count()
    };
    // A's edge x = 2 inside B's top.
    assert_eq!(
        span(
            a_side,
            b_top,
            Point3::new(2.0, 0.5, 0.0),
            Point3::new(2.0, 1.5, 0.0)
        ),
        1,
        "{records:#?}"
    );
    // B's edge x = 1 inside A's bottom, whole.
    assert_eq!(
        span(
            a_bottom,
            b_left,
            Point3::new(1.0, 0.5, 0.0),
            Point3::new(1.0, 1.5, 0.0)
        ),
        1
    );
    // B's edge y = 0.5 inside A's bottom up to A's boundary x = 2.
    assert_eq!(
        span(
            a_bottom,
            b_front,
            Point3::new(1.0, 0.5, 0.0),
            Point3::new(2.0, 0.5, 0.0)
        ),
        1
    );
    assert!(
        !records
            .iter()
            .any(|r| r.face_a == a_bottom && r.face_b == b_top),
        "the coplanar pair is the coincidence engine's"
    );
    // Every record lies in both faces' planes at z = 0.
    for r in &records {
        let (s, e) = ends(r);
        assert!(s.z.abs() < 1e-12 && e.z.abs() < 1e-12, "{r:?}");
    }

    // Records the solver already delivered are not duplicated.
    let n = records.len();
    recover(&a, &b, &mut records, &tol);
    assert_eq!(records.len(), n);
}

/// Edges that run along the receiving face's own boundary divide
/// nothing: two blocks side by side share the line x = 1, z = 0, and
/// neither block's edge on it is imprinted into the other's face.
#[test]
fn an_edge_along_the_faces_own_boundary_is_not_imprinted() {
    let tol = Tolerance::default();
    let a = block((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), tol);
    let b = block((1.0, 0.0, -1.0), (2.0, 1.0, 0.0), tol);
    let mut records = Vec::new();
    recover(&a, &b, &mut records, &tol);
    assert!(records.is_empty(), "{records:#?}");
}

/// A body whose vertices were moved without its curves (a test-built
/// body) has curves where the body is not: they are not read.
#[test]
fn an_edge_whose_curve_misses_its_vertices_is_not_read() {
    let tol = Tolerance::default();
    let a = block((0.0, 0.0, 0.0), (2.0, 2.0, 1.0), tol);
    let mut b = block((1.0, 0.5, -1.0), (3.0, 1.5, 0.0), tol);
    for (_, v) in &mut b.vertices {
        v.position = v.position + Vec3::new(10.0, 0.0, 0.0);
    }
    let mut records = Vec::new();
    recover(&a, &b, &mut records, &tol);
    assert!(
        !records
            .iter()
            .any(|r| (ends(r).0.x - 1.0).abs() < 1e-9 && (ends(r).1.x - 1.0).abs() < 1e-9),
        "B's edge x = 1 is read from a curve its vertices left"
    );
}

/// A wedge whose ridge lies in A's bottom face: either slanted face's
/// record along the ridge imprints the bottom there, so the ridge gets
/// one record, not one per face (the clip gives a cut riding an edge
/// to one of its faces; two copies would have both claim it).
#[test]
fn a_record_of_either_face_along_the_edge_covers_it() {
    let tol = Tolerance::default();
    let a = block((0.0, 0.0, 0.0), (2.0, 2.0, 1.0), tol);
    let profile = [
        Point3::new(0.5, 0.5, -1.0),
        Point3::new(1.5, 0.5, -1.0),
        Point3::new(1.0, 0.5, 0.0),
    ];
    let wedge = BRep::prism(&profile, Vec3::new(0.0, 1.0, 0.0), 1.0, tol).expect("wedge");
    let along_ridge = |r: &IntersectionRecord| {
        runs_between(r, Point3::new(1.0, 0.5, 0.0), Point3::new(1.0, 1.5, 0.0))
    };
    let mut records = Vec::new();
    recover(&a, &wedge, &mut records, &tol);
    assert_eq!(
        records.iter().filter(|r| along_ridge(r)).count(),
        1,
        "{records:#?}"
    );
}
