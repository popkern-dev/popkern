//! Degree elevation for NURBS curves (Piegl & Tiller A5.9).

use super::Homogeneous;
use crate::geometry::Point3;

#[allow(dead_code)]
/// Elevate the degree of a 1D NURBS by 1.
/// Returns (new_degree, new_knots, new_control_points, new_weights).
/// The curve shape is mathematically unchanged.
pub(crate) fn elevate_degree_1d(
    degree: usize,
    knots: &[f64],
    control_points: &[Point3],
    weights: &[f64],
) -> (usize, Vec<f64>, Vec<Point3>, Vec<f64>) {
    let p = degree;
    let new_p = p + 1;

    // Step 1: Find distinct interior knots and their multiplicities
    let mut distinct_interior: Vec<(f64, usize)> = Vec::new();
    {
        let mut i = p + 1;
        while i < knots.len() - p - 1 {
            let knot_val = knots[i];
            let mut mult = 1;
            while i + mult < knots.len() - p - 1 && knots[i + mult] == knot_val {
                mult += 1;
            }
            distinct_interior.push((knot_val, mult));
            i += mult;
        }
    }

    // Step 2: Insert knots to decompose into Bezier segments
    // Each interior knot needs multiplicity = p. Insert until that's achieved.
    let mut knots_to_insert: Vec<f64> = Vec::new();
    for &(knot_val, mult) in &distinct_interior {
        for _ in mult..p {
            knots_to_insert.push(knot_val);
        }
    }

    // Use knot insertion to decompose
    let (_bez_knots, bez_pts, bez_wts) = if knots_to_insert.is_empty() {
        (knots.to_vec(), control_points.to_vec(), weights.to_vec())
    } else {
        // Sort the knots to insert
        knots_to_insert.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        super::knot_insertion::insert_knots_1d(
            degree,
            knots,
            control_points,
            weights,
            &knots_to_insert,
        )
    };

    // Convert Bezier control points to homogeneous
    let bez_hom: Vec<Homogeneous> = bez_pts
        .iter()
        .zip(bez_wts.iter())
        .map(|(pt, &w)| Homogeneous::from_point_weight(pt, w))
        .collect();

    // Step 3: Determine number of Bezier segments
    let n_bez_controls = bez_pts.len();
    let n_segments = (n_bez_controls - 1) / p; // each segment has p+1 control points, sharing endpoints

    if n_segments == 0 {
        // Degenerate case: single point or empty
        return (
            new_p,
            knots.to_vec(),
            control_points.to_vec(),
            weights.to_vec(),
        );
    }

    // Step 4: Elevate each Bezier segment
    let mut elevated_segments: Vec<Vec<Homogeneous>> = Vec::with_capacity(n_segments);
    for seg in 0..n_segments {
        let start = seg * p;
        let seg_pts: Vec<Homogeneous> = (0..=p).map(|i| bez_hom[start + i]).collect();

        // Bezier degree elevation: p -> p+1, producing p+2 control points
        let mut q = vec![
            Homogeneous {
                wx: 0.0,
                wy: 0.0,
                wz: 0.0,
                w: 0.0,
            };
            new_p + 1
        ];
        q[0] = seg_pts[0];
        q[new_p] = seg_pts[p];

        for i in 1..new_p {
            let alpha = i as f64 / new_p as f64;
            q[i] = Homogeneous::lerp(&seg_pts[i], &seg_pts[i - 1], alpha);
        }

        elevated_segments.push(q);
    }

    // Step 5: Merge elevated Bezier segments, removing duplicate boundary points
    let mut new_hom: Vec<Homogeneous> = Vec::new();
    for (seg_idx, seg) in elevated_segments.iter().enumerate() {
        if seg_idx == 0 {
            new_hom.extend_from_slice(seg);
        } else {
            // Skip the first point of subsequent segments (it's the same as the last of previous)
            new_hom.extend_from_slice(&seg[1..]);
        }
    }

    // Step 6: Build new knot vector.
    // The merged Bezier segments join at multiplicity new_p (each pair
    // of adjacent segments shares exactly one control point — a C0
    // representation), so every distinct interior knot carries new_p
    // knots here regardless of its ORIGINAL multiplicity. Writing
    // orig_mult + 1 instead (the pre-2026-08-09 defect) desynchronized
    // the knot vector from the merged control net whenever an interior
    // knot had multiplicity < p (i.e. whenever Step 2 had to insert).
    let mut new_knots: Vec<f64> = Vec::new();

    // Start knot: multiplicity was p+1, becomes p+2 = new_p+1
    let start_knot = knots[0];
    for _ in 0..new_p + 1 {
        new_knots.push(start_knot);
    }

    // Interior knots: segment-join multiplicity (see above)
    for &(knot_val, _) in &distinct_interior {
        for _ in 0..new_p {
            new_knots.push(knot_val);
        }
    }

    // End knot: multiplicity was p+1, becomes p+2 = new_p+1
    let end_knot = knots[knots.len() - 1];
    for _ in 0..new_p + 1 {
        new_knots.push(end_knot);
    }

    // Convert back from homogeneous
    let mut new_pts = Vec::with_capacity(new_hom.len());
    let mut new_wts = Vec::with_capacity(new_hom.len());
    for h in &new_hom {
        let (pt, w) = h.to_point3_and_weight();
        new_pts.push(pt);
        new_wts.push(w);
    }

    // Verify consistency: n_knots = n_controls + degree + 1
    debug_assert_eq!(
        new_knots.len(),
        new_pts.len() + new_p + 1,
        "Knot vector length mismatch: {} knots, {} controls, degree {}",
        new_knots.len(),
        new_pts.len(),
        new_p
    );

    // Step 7: Restore the continuity-economical representation. The
    // elevated curve should carry each interior knot at orig_mult + 1;
    // the segment-join surplus (new_p - orig_mult - 1 knots) reverses
    // an exact insertion, so removal is mathematically exact. Removal
    // is tolerance-gated per knot; if numerical amplification declines
    // it (nano knot spans), the C0-redundant representation is kept —
    // still shape-exact, merely less economical.
    let mut cur = (new_knots, new_pts, new_wts);
    for &(knot_val, orig_mult) in &distinct_interior {
        let surplus = new_p.saturating_sub(orig_mult + 1);
        if surplus == 0 {
            continue;
        }
        let Some(index) = cur
            .0
            .iter()
            .position(|&knot| (knot - knot_val).abs() < 1e-15)
        else {
            continue;
        };
        let (nk, np, nw, _removed) = super::knot_removal::remove_knot_times_1d(
            new_p,
            &cur.0,
            &cur.1,
            &cur.2,
            index,
            surplus,
            ELEVATION_KNOT_RESTORE_TOLERANCE,
        );
        cur = (nk, np, nw);
    }

    (new_p, cur.0, cur.1, cur.2)
}

/// Tolerance for reversing the exact segment-join knot insertions
/// after degree elevation. The removal is mathematically exact; this
/// bound only absorbs floating-point noise in the reconstruction, so
/// it sits well below geometric tolerances.
const ELEVATION_KNOT_RESTORE_TOLERANCE: f64 = 1.0e-10;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Point3;
    use crate::geometry::curve::CurveEval;
    use crate::geometry::nurbs_curve::NurbsCurve;

    #[test]
    fn elevate_degree_1_to_2() {
        // Degree 1 line: should become degree 2 with 3 control points
        let pts = vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)];
        let wts = vec![1.0, 1.0];
        let knots = vec![0.0, 0.0, 1.0, 1.0];

        let (new_deg, new_knots, new_pts, _new_wts) = elevate_degree_1d(1, &knots, &pts, &wts);
        assert_eq!(new_deg, 2);
        assert_eq!(new_pts.len(), 3);
        assert_eq!(new_knots.len(), new_pts.len() + new_deg + 1);
    }

    #[test]
    fn elevate_preserves_line_shape() {
        let pts = vec![Point3::new(0.0, 0.0, 0.0), Point3::new(2.0, 4.0, 0.0)];
        let wts = vec![1.0, 1.0];
        let knots = vec![0.0, 0.0, 1.0, 1.0];

        let (new_deg, new_knots, new_pts, new_wts) = elevate_degree_1d(1, &knots, &pts, &wts);

        let original = NurbsCurve {
            degree: 1,
            knots,
            control_points: pts,
            weights: wts,
        };
        let elevated = NurbsCurve {
            degree: new_deg,
            knots: new_knots,
            control_points: new_pts,
            weights: new_wts,
        };

        for t in [0.0, 0.25, 0.5, 0.75, 0.99] {
            let p_old = original.point_at(t);
            let p_new = elevated.point_at(t);
            assert!((p_old.x - p_new.x).abs() < 1e-10, "Shape changed at t={t}");
            assert!((p_old.y - p_new.y).abs() < 1e-10, "Shape changed at t={t}");
        }
    }

    #[test]
    fn elevate_quadratic_bezier() {
        let pts = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ];
        let wts = vec![1.0, 1.0, 1.0];
        let knots = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];

        let (new_deg, new_knots, new_pts, new_wts) = elevate_degree_1d(2, &knots, &pts, &wts);
        assert_eq!(new_deg, 3);

        let original = NurbsCurve {
            degree: 2,
            knots,
            control_points: pts,
            weights: wts,
        };
        let elevated = NurbsCurve {
            degree: new_deg,
            knots: new_knots,
            control_points: new_pts,
            weights: new_wts,
        };

        for t in [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99] {
            let p_old = original.point_at(t);
            let p_new = elevated.point_at(t);
            assert!((p_old.x - p_new.x).abs() < 1e-10, "x changed at t={t}");
            assert!((p_old.y - p_new.y).abs() < 1e-10, "y changed at t={t}");
        }
    }

    #[test]
    fn elevate_rational_quarter_circle() {
        let w = 1.0_f64 / 2.0_f64.sqrt();
        let pts = vec![
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ];
        let wts = vec![1.0, w, 1.0];
        let knots = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];

        let (new_deg, new_knots, new_pts, new_wts) = elevate_degree_1d(2, &knots, &pts, &wts);

        let _original = NurbsCurve {
            degree: 2,
            knots,
            control_points: pts,
            weights: wts,
        };
        let elevated = NurbsCurve {
            degree: new_deg,
            knots: new_knots,
            control_points: new_pts,
            weights: new_wts,
        };

        // All points should still be on the unit circle
        for t in [0.0, 0.25, 0.5, 0.75, 0.99] {
            let p_new = elevated.point_at(t);
            let r = (p_new.x * p_new.x + p_new.y * p_new.y).sqrt();
            assert!(
                (r - 1.0).abs() < 1e-10,
                "Point at t={t} not on unit circle: r={r}"
            );
        }
    }

    #[test]
    fn elevate_multi_span_rational() {
        // Multi-span rational curve with interior knot of multiplicity 2
        // (already at full Bezier decomposition)
        let w = 1.0_f64 / 2.0_f64.sqrt();
        let pts = vec![
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(-1.0, 1.0, 0.0),
            Point3::new(-1.0, 0.0, 0.0),
        ];
        let wts = vec![1.0, w, 1.0, w, 1.0];
        let knots = vec![0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0];

        let (new_deg, new_knots, new_pts, _new_wts) = elevate_degree_1d(2, &knots, &pts, &wts);
        assert_eq!(new_deg, 3);
        assert_eq!(new_knots.len(), new_pts.len() + new_deg + 1);
    }

    #[test]
    fn elevate_multi_span_simple_interior_knot() {
        // Interior knot BELOW full multiplicity (mult 1 at degree 2):
        // Step 2 must insert to decompose, and the knot vector must
        // stay consistent with the merged control net. This was the
        // pre-2026-08-09 defect: orig_mult + 1 knots were written for
        // a segment join carrying new_p, producing an invalid curve
        // (PK-255: a partially-removed terminal sliver knot fed this
        // path and the fillet exact-boundary preparation panicked).
        let pts = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 2.0, 0.0),
            Point3::new(2.0, 2.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
        ];
        let wts = vec![1.0, 1.0, 1.0, 1.0];
        let knots = vec![0.0, 0.0, 0.0, 0.4, 1.0, 1.0, 1.0];

        let (new_deg, new_knots, new_pts, new_wts) = elevate_degree_1d(2, &knots, &pts, &wts);
        assert_eq!(new_deg, 3);
        assert_eq!(new_knots.len(), new_pts.len() + new_deg + 1);

        let original = NurbsCurve {
            degree: 2,
            knots,
            control_points: pts,
            weights: wts,
        };
        let elevated = NurbsCurve {
            degree: new_deg,
            knots: new_knots,
            control_points: new_pts,
            weights: new_wts,
        };
        for t in [0.0, 0.1, 0.3, 0.4, 0.5, 0.7, 0.9, 1.0] {
            let p_old = original.point_at(t);
            let p_new = elevated.point_at(t);
            assert!(
                p_old.distance_to(&p_new) < 1e-10,
                "Shape changed at t={t}: {p_old:?} vs {p_new:?}"
            );
        }
    }

    #[test]
    fn elevate_near_end_sliver_interior_knot() {
        // The measured PK-255 shape: a simple interior knot 3e-8 from
        // the domain end. Elevation must stay valid and shape-true;
        // if the exact restore removal declines on the nano span, the
        // C0-redundant representation is acceptable.
        let sliver = 1.0 - 2.98e-8;
        let pts = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(2.0, 1.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
        ];
        let wts = vec![1.0, 1.0, 1.0, 1.0];
        let knots = vec![0.0, 0.0, 0.0, sliver, 1.0, 1.0, 1.0];

        let (new_deg, new_knots, new_pts, new_wts) = elevate_degree_1d(2, &knots, &pts, &wts);
        assert_eq!(new_deg, 3);
        assert_eq!(new_knots.len(), new_pts.len() + new_deg + 1);

        let original = NurbsCurve {
            degree: 2,
            knots,
            control_points: pts,
            weights: wts,
        };
        let elevated = NurbsCurve {
            degree: new_deg,
            knots: new_knots,
            control_points: new_pts,
            weights: new_wts,
        };
        for t in [0.0, 0.25, 0.5, 0.75, 0.9999, 1.0] {
            let p_old = original.point_at(t);
            let p_new = elevated.point_at(t);
            assert!(
                p_old.distance_to(&p_new) < 1e-9,
                "Shape changed at t={t}: {p_old:?} vs {p_new:?}"
            );
        }
    }
}
