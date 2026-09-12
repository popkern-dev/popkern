//! Pins for the trim-residue post-filters (PK-487): cross-pair
//! duplicate residues of one span on one owner face collapse to a
//! single copy; distinct spans and distinct owner faces survive.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve_on_surface::CurveOnSurface;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::FaceHandle;

use super::dedupe_same_face_residues;
use crate::intersect::{EndpointBoundaryType, IntersectionRecord};

fn line(a: Point3, b: Point3) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![a, b],
        weights: vec![1.0, 1.0],
    }
}

fn residue(
    face_a: FaceHandle,
    face_b: FaceHandle,
    owner_is_a: bool,
    a: Point3,
    b: Point3,
) -> IntersectionRecord {
    let uv = line(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0));
    IntersectionRecord {
        face_a,
        face_b,
        curve: CurveOnSurface {
            curve_3d: line(a, b),
            curve_uv1: uv.clone(),
            curve_uv2: uv,
            is_closed: false,
        },
        is_tangent: false,
        is_coincident: false,
        start_boundary: EndpointBoundaryType::OnEdgeA,
        end_boundary: EndpointBoundaryType::OnEdgeA,
        residue_owner_a: Some(owner_is_a),
    }
}

fn fh(n: u64) -> FaceHandle {
    // Distinct handles minted through a scratch BRep-free slotmap:
    // FaceHandle is a slotmap key, so build them from a real map.
    let mut map: slotmap::SlotMap<FaceHandle, ()> = slotmap::SlotMap::with_key();
    let mut last = map.insert(());
    for _ in 1..n {
        last = map.insert(());
    }
    last
}

#[test]
fn same_span_same_owner_face_residues_collapse_to_one() {
    let tol = Tolerance::default();
    let owner = fh(1);
    let (p, q) = (Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0));
    // Three pairs each emitted the same span on the same owner-A face,
    // with independent-trace slop on one copy's endpoint.
    let slop = tol.linear() * 20.0;
    let mut records = vec![
        residue(owner, fh(2), true, p, q),
        residue(owner, fh(3), true, p, q),
        residue(owner, fh(4), true, Point3::new(slop, 0.0, 0.0), q),
    ];
    dedupe_same_face_residues(&mut records, &tol);
    assert_eq!(
        records.len(),
        1,
        "same-face same-span residues must collapse to one copy"
    );
}

#[test]
fn distinct_spans_and_distinct_owner_faces_survive() {
    let tol = Tolerance::default();
    let (p, q) = (Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0));
    let r = Point3::new(0.0, 1.0, 0.0);
    let s = Point3::new(1.0, 1.0, 0.0);
    let mut records = vec![
        // Same owner face, DIFFERENT spans: both partition the face.
        residue(fh(1), fh(2), true, p, q),
        residue(fh(1), fh(3), true, r, s),
        // Same span as the first, but a DIFFERENT owner face: it
        // partitions that other face and must survive.
        residue(fh(5), fh(2), false, p, q),
    ];
    dedupe_same_face_residues(&mut records, &tol);
    assert_eq!(
        records.len(),
        3,
        "distinct spans and distinct owner faces are not duplicates"
    );
}
