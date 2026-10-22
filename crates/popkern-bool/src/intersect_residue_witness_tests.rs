//! Pins for the residue closure witness (PK-473): a residue swarm
//! inside an outline the ordinary records close is retired whole, and
//! a chain of more than 16 residues closing an outline the ordinary
//! records leave open survives — the swarm cap it replaces dropped
//! both.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve_on_surface::CurveOnSurface;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::FaceHandle;

use super::{ResidueWitness, retire_unwitnessed_residues};
use crate::intersect::{EndpointBoundaryType, IntersectionRecord};

fn line(a: Point3, b: Point3) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![a, b],
        weights: vec![1.0, 1.0],
    }
}

fn record(
    face_a: FaceHandle,
    face_b: FaceHandle,
    owner: Option<bool>,
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
        residue_owner_a: owner,
    }
}

fn faces(n: usize) -> Vec<FaceHandle> {
    let mut map: slotmap::SlotMap<FaceHandle, ()> = slotmap::SlotMap::with_key();
    (0..n).map(|_| map.insert(())).collect()
}

fn p(x: f64, y: f64) -> Point3 {
    Point3::new(x, y, 0.0)
}

fn gate() -> f64 {
    Tolerance::default().linear() * 20.0
}

/// The unit square's sides as ordinary records on `owner` (face A),
/// each against its own B face; `open_top` leaves the y = 1 side out.
fn square(owner: FaceHandle, others: &[FaceHandle], open_top: bool) -> Vec<IntersectionRecord> {
    let corners = [p(0.0, 0.0), p(1.0, 0.0), p(1.0, 1.0), p(0.0, 1.0)];
    (0..4)
        .filter(|&k| !(open_top && k == 2))
        .map(|k| record(owner, others[k], None, corners[k], corners[(k + 1) % 4]))
        .collect()
}

fn residue_count(records: &[IntersectionRecord]) -> usize {
    records
        .iter()
        .filter(|r| r.residue_owner_a.is_some())
        .count()
}

#[test]
fn a_chord_swarm_inside_a_closed_outline_retires_whole() {
    let f = faces(6);
    let mut records = square(f[0], &f[1..5], false);
    // 24 crisscrossing chords across the square, the impeller floor's
    // shape: every end on an outline side's interior.
    for k in 0..24 {
        let s = (k as f64 + 0.5) / 24.0;
        records.push(record(f[0], f[5], Some(true), p(s, 0.0), p(1.0 - s, 1.0)));
    }
    let rows = retire_unwitnessed_residues(&mut records, gate(), |_, _, _| false);
    assert_eq!(residue_count(&records), 0);
    assert_eq!(records.len(), 4, "ordinary records are never retired");
    assert_eq!(
        rows,
        vec![ResidueWitness {
            owner_is_a: true,
            face: f[0],
            dangling_ends: 0,
            kept: 0,
            retired: 24
        }]
    );
}

#[test]
fn a_chain_of_twenty_residues_closing_an_open_outline_survives() {
    let f = faces(6);
    // The y = 1 side is missing (a graze the intersector lost): the
    // outline's ends at (1, 1) and (0, 1) dangle.
    let mut records = square(f[0], &f[1..5], true);
    for k in 0..20 {
        let (x0, x1) = (1.0 - k as f64 / 20.0, 1.0 - (k + 1) as f64 / 20.0);
        records.push(record(f[0], f[5], Some(true), p(x0, 1.0), p(x1, 1.0)));
    }
    // Plus chords that only subdivide: retired.
    for k in 0..5 {
        let s = (k as f64 + 0.5) / 5.0;
        records.push(record(f[0], f[5], Some(true), p(0.0, s), p(1.0, s)));
    }
    let rows = retire_unwitnessed_residues(&mut records, gate(), |_, _, _| false);
    assert_eq!(rows.len(), 1);
    assert_eq!(
        (rows[0].dangling_ends, rows[0].kept, rows[0].retired),
        (2, 20, 5)
    );
    // The survivors are exactly the y = 1 chain.
    assert!(
        records
            .iter()
            .filter(|r| r.residue_owner_a.is_some())
            .all(|r| r.curve.curve_3d.control_points.iter().all(|c| c.y == 1.0))
    );
}

#[test]
fn ordinary_ends_on_the_owner_boundary_do_not_dangle() {
    let f = faces(3);
    // One ordinary record crossing the face from boundary to boundary
    // (x = 0 and x = 1 are the face's trim), with a residue continuing
    // from each end along the boundary.
    let mut records = vec![
        record(f[0], f[1], None, p(0.0, 0.5), p(1.0, 0.5)),
        record(f[0], f[2], Some(true), p(1.0, 0.5), p(1.0, 1.0)),
        record(f[0], f[2], Some(true), p(0.0, 0.5), p(0.0, 0.0)),
    ];
    let on_boundary = |_: bool, _: FaceHandle, q: &Point3| q.x == 0.0 || q.x == 1.0;
    let rows = retire_unwitnessed_residues(&mut records, gate(), on_boundary);
    assert_eq!((rows[0].dangling_ends, rows[0].retired), (0, 2));
    // Without the boundary, both ends dangle and both residues stay.
    let mut records = vec![
        record(f[0], f[1], None, p(0.0, 0.5), p(1.0, 0.5)),
        record(f[0], f[2], Some(true), p(1.0, 0.5), p(1.0, 1.0)),
        record(f[0], f[2], Some(true), p(0.0, 0.5), p(0.0, 0.0)),
    ];
    let rows = retire_unwitnessed_residues(&mut records, gate(), |_, _, _| false);
    assert_eq!((rows[0].dangling_ends, rows[0].kept), (2, 2));
}

#[test]
fn a_dangling_end_on_one_face_witnesses_no_other_face() {
    let f = faces(4);
    // An open outline on A face f0; a residue owned by B face f1 runs
    // from the same dangling end. It partitions f1, whose ordinary
    // records (none here) leave nothing open: retired. The A residue
    // at the same spot survives.
    let mut records = square(f[0], &[f[1], f[2], f[3], f[2]], true);
    records.push(record(f[0], f[1], Some(false), p(1.0, 1.0), p(0.0, 1.0)));
    records.push(record(f[0], f[3], Some(true), p(1.0, 1.0), p(0.0, 1.0)));
    let rows = retire_unwitnessed_residues(&mut records, gate(), |_, _, _| false);
    let row = |a: bool| {
        rows.iter()
            .find(|r| r.owner_is_a == a)
            .map(|r| (r.kept, r.retired))
    };
    assert_eq!(row(false), Some((0, 1)));
    assert_eq!(row(true), Some((1, 0)));
}
