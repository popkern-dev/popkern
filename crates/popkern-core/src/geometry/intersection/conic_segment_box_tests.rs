//! PK-785: conic segment boxes contain the segment and are tight.

use super::exact_conic_segment_box;
use crate::geometry::curve::{CurveEval, CurveGeometry};
use crate::geometry::{Point3, Vec3};

fn conics() -> Vec<CurveGeometry> {
    let normal = Vec3::new(0.3, -0.5, 0.81)
        .normalized()
        .expect("unit normal");
    let center = Point3::new(1.0, -2.0, 0.5);
    vec![
        CurveGeometry::Circle {
            center,
            normal,
            radius: 5.0,
        },
        CurveGeometry::Arc {
            center,
            normal,
            radius: 5.0,
            start_angle: -0.5,
            end_angle: 0.5,
        },
        CurveGeometry::Arc {
            center,
            normal,
            radius: 2.0,
            start_angle: 5.9,
            end_angle: 7.4,
        },
        CurveGeometry::Ellipse {
            center,
            u_basis: Vec3::new(3.0, 0.0, 0.0),
            v_basis: Vec3::new(0.0, 1.2, 0.9),
            start_angle: -2.0,
            end_angle: 2.5,
        },
    ]
}

/// Every Bézier segment of each conic: dense samples inside the box,
/// and each box face reached by the samples (no over-cover beyond the
/// sample spacing's sagitta).
#[test]
fn conic_segment_boxes_are_exact() {
    for conic in conics() {
        let (t0, t1) = match &conic {
            CurveGeometry::Circle { .. } => (0.0, std::f64::consts::TAU),
            CurveGeometry::Arc {
                start_angle,
                end_angle,
                ..
            }
            | CurveGeometry::Ellipse {
                start_angle,
                end_angle,
                ..
            } => (*start_angle, *end_angle),
            _ => unreachable!(),
        };
        let nurbs = conic.to_nurbs_in_range((t0, t1));
        let segments = super::super::super::decompose_to_bezier_segments(&nurbs);
        assert!(!segments.is_empty());
        for (seg, tr) in segments {
            let Some((lo, hi)) = exact_conic_segment_box(&conic, &seg, tr) else {
                panic!("{conic:?} segment {tr:?} has an exact box");
            };
            let n = 2000;
            let mut seen_lo = [f64::INFINITY; 3];
            let mut seen_hi = [f64::NEG_INFINITY; 3];
            for i in 0..=n {
                let p = seg.point_at(tr.0 + (tr.1 - tr.0) * i as f64 / n as f64);
                for (k, c) in [p.x, p.y, p.z].into_iter().enumerate() {
                    let (l, h) = ([lo.x, lo.y, lo.z][k], [hi.x, hi.y, hi.z][k]);
                    assert!(
                        c >= l && c <= h,
                        "{conic:?} {tr:?} axis {k}: {c} outside [{l}, {h}]"
                    );
                    seen_lo[k] = seen_lo[k].min(c);
                    seen_hi[k] = seen_hi[k].max(c);
                }
            }
            for k in 0..3 {
                let (l, h) = ([lo.x, lo.y, lo.z][k], [hi.x, hi.y, hi.z][k]);
                assert!(
                    seen_lo[k] - l < 1e-5 && h - seen_hi[k] < 1e-5,
                    "{conic:?} {tr:?} axis {k}: box [{l}, {h}] vs samples [{}, {}]",
                    seen_lo[k],
                    seen_hi[k]
                );
            }
        }
    }
}

#[test]
fn other_curves_keep_their_sampled_boxes() {
    let line = CurveGeometry::Line {
        origin: Point3::new(0.0, 0.0, 0.0),
        direction: Vec3::new(1.0, 0.0, 0.0),
    };
    let nurbs = line.to_nurbs_in_range((0.0, 1.0));
    assert!(exact_conic_segment_box(&line, &nurbs, (0.0, 1.0)).is_none());
}
