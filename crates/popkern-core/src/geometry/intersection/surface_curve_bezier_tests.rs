//! One-pass Bezier decomposition agrees with knot-by-knot splitting.

use super::{
    bezier_segments_by_refinement, bezier_segments_by_splitting, decompose_to_bezier_segments,
    interior_knots_and_insertions,
};
use crate::geometry::Point3;
use crate::geometry::nurbs_curve::NurbsCurve;

/// A clamped rational curve with `spans` spans, interior knot
/// multiplicities cycling 1..=degree.
fn curve(degree: usize, spans: usize, seed: u64) -> NurbsCurve {
    let mut state = seed;
    let mut unit = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    let mut knots = vec![0.0; degree + 1];
    for s in 1..spans {
        let mult = 1 + s % degree;
        knots.extend(std::iter::repeat_n(s as f64 + 0.3 * unit(), mult));
    }
    knots.extend(std::iter::repeat_n(spans as f64, degree + 1));
    let count = knots.len() - degree - 1;
    NurbsCurve {
        degree,
        knots,
        control_points: (0..count)
            .map(|i| Point3::new(i as f64 + unit(), 3.0 * unit(), unit() - 0.5))
            .collect(),
        weights: (0..count).map(|_| 0.5 + unit()).collect(),
    }
}

#[test]
fn refinement_slices_match_knot_by_knot_splitting() {
    for (degree, spans, seed) in [
        (1, 6, 1),
        (2, 9, 2),
        (3, 14, 3),
        (4, 7, 4),
        (5, 11, 5),
        (3, 40, 6),
    ] {
        let c = curve(degree, spans, seed);
        let (interior, inserts) = interior_knots_and_insertions(&c);
        assert!(
            bezier_segments_by_refinement(&c, &interior, &inserts).is_some(),
            "degree {degree}: the refinement path must serve a clamped curve"
        );
        let fast = decompose_to_bezier_segments(&c);
        let reference = bezier_segments_by_splitting(&c);
        assert_eq!(fast.len(), spans, "degree {degree}");
        assert_eq!(fast.len(), reference.len(), "degree {degree}");
        for ((a, ra), (b, rb)) in fast.iter().zip(&reference) {
            assert_eq!(ra, rb, "degree {degree}: segment range");
            assert_eq!(a.degree, b.degree);
            assert_eq!(a.knots, b.knots, "degree {degree}: segment knots");
            assert_eq!(a.control_points.len(), b.control_points.len());
            for (p, q) in a.control_points.iter().zip(&b.control_points) {
                assert!(p.distance_to(q) < 1e-12, "degree {degree}: {p:?} vs {q:?}");
            }
            for (w, v) in a.weights.iter().zip(&b.weights) {
                assert!((w - v).abs() < 1e-12, "degree {degree}: weight {w} vs {v}");
            }
        }
    }
}

/// Knots 1e-15 apart are one knot to the decomposition but two to the
/// refinement: the refined vector is not Bezier form, so the reference
/// path serves and the output is the reference's.
#[test]
fn a_near_duplicate_knot_takes_the_reference_path() {
    let mut c = curve(3, 5, 7);
    let k = c.knots[4];
    c.knots.insert(5, k + 1e-15);
    c.control_points.push(Point3::new(9.0, 0.0, 0.0));
    c.weights.push(1.0);
    let (interior, inserts) = interior_knots_and_insertions(&c);
    assert!(bezier_segments_by_refinement(&c, &interior, &inserts).is_none());
    let fast = decompose_to_bezier_segments(&c);
    let reference = bezier_segments_by_splitting(&c);
    assert_eq!(fast.len(), reference.len());
    for ((a, ra), (b, rb)) in fast.iter().zip(&reference) {
        assert_eq!(ra, rb);
        assert_eq!(a.knots, b.knots);
        assert_eq!(a.control_points, b.control_points);
    }
}
