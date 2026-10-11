//! One-pass knot refinement equals repeated single insertion.

use super::*;
use crate::geometry::nurbs_ops::knot_insertion::insert_knot_1d;

/// A rational cubic with interior knots of multiplicity one and two.
fn curve() -> (usize, Vec<f64>, Vec<Point3>, Vec<f64>) {
    let knots = vec![0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 3.5, 4.0, 4.0, 4.0, 4.0];
    let points: Vec<Point3> = (0..8)
        .map(|i| {
            let t = i as f64;
            Point3::new(t, (0.7 * t).sin() * 3.0, (t * t) * 0.1 - 1.0)
        })
        .collect();
    let weights = vec![1.0, 0.8, 1.3, 1.0, 0.6, 1.1, 0.9, 1.0];
    (3, knots, points, weights)
}

#[test]
fn refinement_matches_repeated_insertion() {
    let (p, knots, points, weights) = curve();
    for inserts in [
        vec![0.5, 1.5, 2.7, 3.75],
        vec![0.0, 0.25, 0.25, 2.0, 4.0],
        vec![3.9],
        vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
    ] {
        let (rk, rp, rw) =
            refine_knot_vector_points(p, &knots, &points, &weights, &inserts).expect("refines");
        let (mut k, mut pts, mut wts) = (knots.clone(), points.clone(), weights.clone());
        for &x in &inserts {
            (k, pts, wts) = insert_knot_1d(p, &k, &pts, &wts, x);
        }
        assert_eq!(rk, k, "knots for {inserts:?}");
        assert_eq!(rp.len(), pts.len());
        for (a, b) in rp.iter().zip(&pts) {
            assert!(a.distance_to(b) < 1e-12, "{inserts:?}: {a:?} vs {b:?}");
        }
        for (a, b) in rw.iter().zip(&wts) {
            assert!((a - b).abs() < 1e-12, "{inserts:?}: weight {a} vs {b}");
        }
    }
}

#[test]
fn invalid_inputs_are_refused() {
    let (p, knots, points, weights) = curve();
    assert!(refine_knot_vector_points(p, &knots, &points, &weights, &[4.5]).is_none());
    assert!(refine_knot_vector_points(p, &knots, &points, &weights, &[2.0, 1.0]).is_none());
    assert!(refine_knot_vector_points(p, &knots[1..], &points, &weights, &[1.0]).is_none());
    let (k, pts, _) =
        refine_knot_vector_points(p, &knots, &points, &weights, &[]).expect("no inserts");
    assert_eq!(k, knots);
    for (a, b) in pts.iter().zip(&points) {
        assert!(a.distance_to(b) < 1e-12, "{a:?} vs {b:?}");
    }
}
