//! The record-order rule (plan §4.3a gate): a vertex never reverses the
//! order of points along a record.

use popkern_core::geometry::Point3;
use popkern_core::topology::VertexHandle;

use super::candidates::{CandidateVertex, Incidence, SourceKind};
use super::resolve::{ResolutionContext, VertexGroup, resolve};
use super::tests::{TOL, cube, line, options, record, vertex};

/// `c` lying on record 0 at `s` too.
fn on_record_0(mut c: CandidateVertex, s: f64) -> CandidateVertex {
    c.incidence.push(Incidence::Record {
        record: 0,
        s,
        end: None,
    });
    c
}

/// A point of record 0's interior at `s`.
fn interior(s: f64, at: Point3) -> CandidateVertex {
    CandidateVertex {
        position: at,
        radius: TOL,
        zone: TOL,
        source: SourceKind::EdgeFaceHit,
        incidence: vec![Incidence::Record {
            record: 0,
            s,
            end: None,
        }],
        key: (SourceKind::EdgeFaceHit, 0, 0, 0, 0),
    }
}

/// Two vertices of one operand 0.6 tol apart on a record, and a third
/// point 0.6 tol beyond the second, within reach of both: joining the
/// first vertex would run the record V → W → V, so it joins the second.
/// The partition offers a candidate the vertex nearest it first
/// (PK-1165), so the point joins the second without trying the first;
/// the order rule itself is the next test's.
#[test]
fn a_point_never_joins_a_vertex_across_another_along_a_record() {
    let c = cube();
    let vs: Vec<VertexHandle> = c.vertices.keys().take(2).collect();
    let records = [record(line(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
    ))];
    let x = |k: f64| 1.0 + k * TOL;
    let cands = vec![
        on_record_0(vertex(0, vs[0], Point3::new(x(0.0), 0.0, 0.0)), 0.5),
        on_record_0(
            vertex(0, vs[1], Point3::new(x(0.6), 0.0, 0.0)),
            0.5 * x(0.6),
        ),
        interior(0.5 * x(1.2), Point3::new(x(1.2), 0.0, 0.0)),
    ];
    let ctx = ResolutionContext {
        operands: [&c, &c],
        records: &records,
    };
    let res = resolve(&cands, &ctx, options(10.0 * TOL)).expect("resolves");
    assert_eq!(res.groups.len(), 2, "{res:?}");
    assert_eq!(res.group_of[2], res.group_of[1], "{res:?}");
    assert!(res.order_conflicts.is_empty(), "{:?}", res.order_conflicts);
}

/// A point between two members of one vertex along a record lies inside
/// its span; a point beyond them, or on a closed record, does not.
#[test]
fn a_point_between_a_vertexs_members_along_a_record_is_inside_its_span() {
    let open = [record(line(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
    ))];
    let cands = vec![
        interior(0.50, Point3::new(1.0, 0.0, 0.0)),
        interior(0.52, Point3::new(1.04, 0.0, 0.0)),
        interior(0.51, Point3::new(1.02, 0.0, 0.0)),
        interior(0.53, Point3::new(1.06, 0.0, 0.0)),
    ];
    let groups = [VertexGroup {
        anchor: 0,
        members: vec![0, 1],
        tolerance: TOL,
    }];
    let inside =
        |records: &[_], c| super::record_order::reversal(&cands, records, &groups, None, c);
    assert_eq!(inside(&open, 2), Some((0, 0)));
    assert_eq!(inside(&open, 3), None);
    // Joining a vertex with no point on the record does not take a
    // point out of the span it lies in.
    let off_record = CandidateVertex {
        incidence: vec![],
        ..interior(0.0, Point3::new(1.02, 0.0, 0.0))
    };
    let mut cands = cands.clone();
    cands.push(off_record);
    let groups = [
        groups[0].clone(),
        VertexGroup {
            anchor: 4,
            members: vec![4],
            tolerance: TOL,
        },
    ];
    assert_eq!(
        super::record_order::reversal(&cands, &open, &groups, Some(1), 2),
        Some((0, 0))
    );
    let mut closed = open.clone();
    closed[0].curve.is_closed = true;
    assert_eq!(inside(&closed, 2), None);
}
