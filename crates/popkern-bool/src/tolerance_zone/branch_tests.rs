//! The branch rule (plan §4.3a gate): two branches of one face pair
//! passing within tolerance stay two vertices; two pieces of one branch
//! meeting end to end, and a record ending on another face pair's
//! record, are one.

use popkern_core::geometry::Point3;
use popkern_core::topology::FaceHandle;

use super::candidates::{CandidateVertex, Incidence, SourceKind};
use super::resolve::{Incompatibility, ResolutionContext, resolve};
use super::tests::{TOL, cube, line, options, record, record_end};
use crate::intersect::IntersectionRecord;

/// Two face handles of one operand, for records of distinct pairs.
fn faces() -> (FaceHandle, FaceHandle) {
    let c = cube();
    let mut keys = c.faces.keys();
    let a = keys.next().expect("face");
    let b = keys.next().expect("face");
    (a, b)
}

fn on_pair(mut r: IntersectionRecord, face_a: FaceHandle) -> IntersectionRecord {
    r.face_a = face_a;
    r
}

/// A point on record `record`'s interior at `s` (an edge-face hit the
/// record passes through, say).
fn on_record(record: usize, s: f64, at: Point3) -> CandidateVertex {
    CandidateVertex {
        position: at,
        radius: TOL,
        zone: TOL,
        source: SourceKind::EdgeFaceHit,
        incidence: vec![Incidence::Record {
            record,
            s,
            end: None,
        }],
        key: (SourceKind::EdgeFaceHit, 0, record as u64, 0, 0),
    }
}

/// A record of one face pair ending 0.5 tol beside another record of
/// that pair, mid-span: two branches passing close. The SSI marched
/// the first through, seeing no branch point, so the end is not on it
/// (no T-junction) and the point it passes stays a vertex of its own.
#[test]
fn two_branches_half_a_tolerance_apart_stay_two_vertices() {
    let c = cube().transform(&popkern_core::geometry::transform::Transform3::translation(
        100.0, 0.0, 0.0,
    ));
    let c = c.expect("cube away");
    let records = [
        record(line(Point3::new(0.0, 0.0, 0.0), Point3::new(2.0, 0.0, 0.0))),
        record(line(
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(1.0, 0.5 * TOL, 0.0),
        )),
    ];
    let mut cands = vec![
        record_end(1, 1.0, true, Point3::new(1.0, 0.5 * TOL, 0.0)),
        on_record(0, 0.5, Point3::new(1.0, 0.0, 0.0)),
    ];
    super::incidence::add_interior_incidences(&mut cands, [&c, &c], &records, TOL);
    assert!(
        !super::incidence::is_t_junction(&cands[0]),
        "a branch's end is not on another branch: {:?}",
        cands[0].incidence
    );
    let ctx = ResolutionContext {
        operands: [&c, &c],
        records: &records,
    };
    let res = resolve(&cands, &ctx, options(10.0 * TOL)).expect("resolves");
    assert_eq!(res.groups.len(), 2, "{res:?}");
    assert!(
        res.separations
            .iter()
            .any(|&(_, _, why)| why == Incompatibility::DistinctBranches),
        "{:?}",
        res.separations
    );
}

/// Two records of one face pair meeting end to end (one branch cut at
/// a seam or a band lattice) are one vertex; so is a record ending
/// mid-span of a record of another face pair (a T-junction).
#[test]
fn pieces_of_one_branch_and_junctions_of_two_pairs_are_one_vertex() {
    let c = cube();
    let (fa, fb) = faces();
    let pieces = [
        record(line(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0))),
        record(line(
            Point3::new(1.0, 0.5 * TOL, 0.0),
            Point3::new(2.0, 1.0, 0.0),
        )),
    ];
    let cands = vec![
        record_end(0, 1.0, true, Point3::new(1.0, 0.0, 0.0)),
        record_end(1, 0.0, false, Point3::new(1.0, 0.5 * TOL, 0.0)),
    ];
    let ctx = ResolutionContext {
        operands: [&c, &c],
        records: &pieces,
    };
    let res = resolve(&cands, &ctx, options(10.0 * TOL)).expect("resolves");
    assert_eq!(res.groups.len(), 1, "pieces of one branch: {res:?}");

    let pairs = [
        on_pair(
            record(line(Point3::new(0.0, 0.0, 0.0), Point3::new(2.0, 0.0, 0.0))),
            fa,
        ),
        on_pair(
            record(line(
                Point3::new(1.0, 1.0, 0.0),
                Point3::new(1.0, 0.5 * TOL, 0.0),
            )),
            fb,
        ),
    ];
    let mut cands = vec![
        record_end(1, 1.0, true, Point3::new(1.0, 0.5 * TOL, 0.0)),
        on_record(0, 0.5, Point3::new(1.0, 0.0, 0.0)),
    ];
    let far = c
        .clone()
        .transform(&popkern_core::geometry::transform::Transform3::translation(
            100.0, 0.0, 0.0,
        ))
        .expect("cube away");
    super::incidence::add_interior_incidences(&mut cands, [&far, &far], &pairs, TOL);
    assert!(super::incidence::is_t_junction(&cands[0]));
    let ctx = ResolutionContext {
        operands: [&far, &far],
        records: &pairs,
    };
    let res = resolve(&cands, &ctx, options(10.0 * TOL)).expect("resolves");
    assert_eq!(res.groups.len(), 1, "a junction of two pairs: {res:?}");
}
