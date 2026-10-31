//! Pave blocks and shrunk ranges (plan §4.3a gate): a piece within its
//! vertices is a micro block; two pieces within tolerance of each other
//! between the same two vertices are a common block.

use popkern_core::geometry::Point3;
use popkern_core::topology::VertexHandle;

use super::pave_blocks::{BlockOwner, common_block_pairs, pave_blocks};
use super::resolve::{ResolutionContext, resolve};
use super::tests::{TOL, cube, line, options, quadratic, record, record_end, vertex};
use crate::intersect::IntersectionRecord;
use popkern_core::geometry::nurbs_curve::NurbsCurve;

/// Two vertices of one operand on a record's interior stay two (a
/// valid operand's vertices are distinct). 0.8 tol apart the piece
/// between them never leaves either ball; 1.5 tol apart it leaves each
/// only inside the other: both are micro blocks. 2.5 tol apart the
/// piece is a curve of its own, as are the pieces out to the ends.
#[test]
fn a_piece_between_two_vertices_within_both_is_a_micro_block() {
    let c = cube();
    let vs: Vec<VertexHandle> = c.vertices.keys().take(2).collect();
    let records = [record(line(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
    ))];
    let on_record = |mut cand: super::candidates::CandidateVertex, s: f64| {
        cand.incidence.push(super::candidates::Incidence::Record {
            record: 0,
            s,
            end: None,
        });
        cand
    };
    for (gap, micro) in [(0.8, true), (1.5, true), (2.5, false)] {
        let x = 0.5 + gap * TOL;
        let cands = vec![
            on_record(vertex(0, vs[0], Point3::new(0.5, 0.0, 0.0)), 0.5),
            on_record(vertex(0, vs[1], Point3::new(x, 0.0, 0.0)), x),
            record_end(0, 0.0, false, Point3::new(0.0, 0.0, 0.0)),
            record_end(0, 1.0, true, Point3::new(1.0, 0.0, 0.0)),
        ];
        let ctx = ResolutionContext {
            operands: [&c, &c],
            records: &records,
        };
        let res = resolve(&cands, &ctx, options(10.0 * TOL)).expect("resolves");
        let blocks = pave_blocks(&cands, &res, [&c, &c], &records, TOL);
        assert_eq!(blocks.len(), 3, "gap {gap}: {blocks:?}");
        assert!(blocks.iter().all(|b| b.owner == BlockOwner::Record(0)));
        let found: Vec<_> = blocks.iter().filter(|b| b.shrunk.is_none()).collect();
        assert_eq!(found.len(), usize::from(micro), "gap {gap}: {blocks:?}");
        if let [m] = found.as_slice() {
            assert_eq!(m.ends, [Some(res.group_of[0]), Some(res.group_of[1])]);
        }
    }
}

fn two_records(
    second: NurbsCurve,
) -> (
    [IntersectionRecord; 2],
    Vec<super::candidates::CandidateVertex>,
) {
    let records = [
        record(line(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0))),
        record(second),
    ];
    let cands = vec![
        record_end(0, 0.0, false, Point3::new(0.0, 0.0, 0.0)),
        record_end(0, 1.0, true, Point3::new(1.0, 0.0, 0.0)),
        record_end(1, 0.0, false, Point3::new(0.0, 0.0, 0.0)),
        record_end(1, 1.0, true, Point3::new(1.0, 0.0, 0.0)),
    ];
    (records, cands)
}

/// Two records between the same two vertices, half a tolerance apart
/// all along, trace one piece: a common block (served twice, a face
/// no wider than tolerance). Bowed 5 tol apart in the middle they bound
/// a face of their own.
#[test]
fn two_pieces_within_tolerance_between_the_same_vertices_are_a_common_block() {
    let c = cube();
    for (bow, common) in [(0.5 * TOL, true), (10.0 * TOL, false)] {
        let (records, cands) = two_records(quadratic([
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, bow, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ]));
        let ctx = ResolutionContext {
            operands: [&c, &c],
            records: &records,
        };
        let res = resolve(&cands, &ctx, options(10.0 * TOL)).expect("resolves");
        assert_eq!(res.groups.len(), 2, "{res:?}");
        let blocks = pave_blocks(&cands, &res, [&c, &c], &records, TOL);
        assert_eq!(blocks.len(), 2, "{blocks:?}");
        let found = common_block_pairs(&blocks, &res, [&c, &c], &records);
        assert_eq!(found.len(), usize::from(common), "bow {bow}: {found:?}");
    }
}
