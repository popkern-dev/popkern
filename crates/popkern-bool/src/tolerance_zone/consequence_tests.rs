//! Consequences of a resolution (plan §4.3a): a record whose two ends
//! resolved into one vertex is a micro edge, a spur or a loop, by the
//! vertex's own tolerance; a record whose ends stay apart, or a closed
//! record, is none of these.

use popkern_core::geometry::Point3;
use popkern_core::geometry::nurbs_curve::NurbsCurve;

use super::consequences::{Collapse, collapsed_ends};
use super::resolve::{ResolutionContext, resolve};
use super::tests::{TOL, cube, line, options, quadratic, record, record_end};

/// A cubic through `points`.
fn cubic(points: [Point3; 4]) -> NurbsCurve {
    NurbsCurve {
        degree: 3,
        knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        control_points: points.to_vec(),
        weights: vec![1.0; 4],
    }
}

/// A record shorter than its ends' tolerances lies within their vertex:
/// a micro edge, which collapses. One running a millimetre out and back
/// along a line is a spur (measured: SSI noise where a wall plane
/// grazes a rim fillet), and one enclosing width is a loop. The
/// curl-back keeps its ends apart and is none, and a closed record's
/// ends meet by construction: no consequence either.
#[test]
fn a_record_whose_ends_meet_is_a_micro_edge_a_spur_or_a_loop() {
    let c = cube();
    let p = |x: f64, y: f64| Point3::new(x, y, 0.0);
    let mut closed = record(quadratic([p(9.0, 0.0), p(11.0, 0.0), p(9.0, 0.0)]));
    closed.curve.is_closed = true;
    let records = [
        record(line(p(0.0, 0.0), p(0.5 * TOL, 0.0))),
        record(quadratic([p(5.0, 0.0), p(7.0, 0.0), p(5.0, 0.5 * TOL)])),
        closed,
        record(cubic([
            p(20.0, 0.0),
            p(20.0, 1e-3),
            p(20.0, -2e-3),
            p(20.0, 0.0),
        ])),
        record(cubic([
            p(30.0, 0.0),
            p(31.0, 1.0),
            p(31.0, -1.0),
            p(30.0, 0.0),
        ])),
    ];
    let ctx = ResolutionContext {
        operands: [&c, &c],
        records: &records,
    };
    let mut cands = Vec::new();
    for (ri, (s, e)) in [
        (p(0.0, 0.0), p(0.5 * TOL, 0.0)),
        (p(5.0, 0.0), p(5.0, 0.5 * TOL)),
        (p(9.0, 0.0), p(9.0, 0.0)),
        (p(20.0, 0.0), p(20.0, 0.0)),
        (p(30.0, 0.0), p(30.0, 0.0)),
    ]
    .into_iter()
    .enumerate()
    {
        cands.push(record_end(ri, 0.0, false, s));
        cands.push(record_end(ri, 1.0, true, e));
    }
    let res = resolve(&cands, &ctx, options(20.0 * TOL)).expect("resolves");
    let kinds: Vec<(usize, Collapse)> = collapsed_ends(&cands, &res, &records)
        .iter()
        .map(|c| (c.record, c.kind()))
        .collect();
    assert_eq!(
        kinds,
        vec![
            (0, Collapse::MicroEdge),
            (3, Collapse::Spur),
            (4, Collapse::Loop)
        ]
    );
}
