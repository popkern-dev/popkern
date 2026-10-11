//! Bezier decomposition of a NURBS curve for the surface x curve grid.

use crate::geometry::curve::CurveEval;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::nurbs_ops::knot_refinement::refine_knot_vector_points;

/// Decompose a NURBS curve into Bezier segments, one per knot span,
/// each with its parameter range.
///
/// Every interior knot is raised to multiplicity `degree` in one knot
/// refinement pass (NURBS Book A5.4), and each segment's control points
/// are then a `degree + 1` slice of the refined polygon. Raising and
/// splitting one knot at a time was quadratic in the span count: on an
/// 80 mm thread it was a third of the clip's endpoint refinement
/// (PK-794). The two agree to rounding (the single insertions
/// round-trip every point through homogeneous form). A curve whose
/// refined knot vector is not exactly clamped Bezier form, such as an
/// unclamped curve or knots closer than 1e-14 but not equal, takes the
/// reference path.
pub(super) fn decompose_to_bezier_segments(nurbs: &NurbsCurve) -> Vec<(NurbsCurve, (f64, f64))> {
    let (interior_knots, knots_to_insert) = interior_knots_and_insertions(nurbs);
    if interior_knots.is_empty() {
        // Already a single Bezier segment.
        return vec![(nurbs.clone(), nurbs.parameter_range())];
    }
    bezier_segments_by_refinement(nurbs, &interior_knots, &knots_to_insert)
        .unwrap_or_else(|| bezier_segments_by_splitting(nurbs))
}

/// The distinct interior knots (values within 1e-14 are one knot) and
/// the insertions that raise each to multiplicity `degree`.
fn interior_knots_and_insertions(nurbs: &NurbsCurve) -> (Vec<f64>, Vec<f64>) {
    let (t_min, t_max) = nurbs.parameter_range();
    let mut interior_knots: Vec<f64> = Vec::new();
    for &k in &nurbs.knots {
        if k > t_min
            && k < t_max
            && interior_knots
                .last()
                .is_none_or(|&last| (k - last).abs() > 1e-14)
        {
            interior_knots.push(k);
        }
    }
    let mut knots_to_insert: Vec<f64> = Vec::new();
    for &k in &interior_knots {
        let mult = nurbs
            .knots
            .iter()
            .filter(|&&v| (v - k).abs() < 1e-14)
            .count();
        for _ in 0..nurbs.degree.saturating_sub(mult) {
            knots_to_insert.push(k);
        }
    }
    (interior_knots, knots_to_insert)
}

/// Segments sliced from the refined polygon, or `None` when the refined
/// knot vector is not `[t_min x (p+1), k_1 x p, ..., k_m x p, t_max x (p+1)]`
/// exactly.
fn bezier_segments_by_refinement(
    nurbs: &NurbsCurve,
    interior_knots: &[f64],
    knots_to_insert: &[f64],
) -> Option<Vec<(NurbsCurve, (f64, f64))>> {
    let p = nurbs.degree;
    if p == 0 {
        return None;
    }
    let (t_min, t_max) = nurbs.parameter_range();
    let (knots, points, weights) = refine_knot_vector_points(
        p,
        &nurbs.knots,
        &nurbs.control_points,
        &nurbs.weights,
        knots_to_insert,
    )?;
    let spans = interior_knots.len() + 1;
    let mut breaks = Vec::with_capacity(spans + 1);
    breaks.push(t_min);
    breaks.extend_from_slice(interior_knots);
    breaks.push(t_max);
    let mut expected = vec![t_min; p + 1];
    for &k in interior_knots {
        expected.extend(std::iter::repeat_n(k, p));
    }
    expected.extend(std::iter::repeat_n(t_max, p + 1));
    if points.len() != spans * p + 1 || knots != expected {
        return None;
    }
    Some(
        breaks
            .windows(2)
            .enumerate()
            .map(|(i, w)| {
                let slice = i * p..=i * p + p;
                let mut segment_knots = vec![w[0]; p + 1];
                segment_knots.extend(std::iter::repeat_n(w[1], p + 1));
                let segment = NurbsCurve {
                    degree: p,
                    knots: segment_knots,
                    control_points: points[slice.clone()].to_vec(),
                    weights: weights[slice].to_vec(),
                };
                (segment, (w[0], w[1]))
            })
            .collect(),
    )
}

/// The reference decomposition: raise every interior knot to multiplicity
/// `degree` by single insertions, then split off one segment at a time.
/// Both steps rebuild the whole remaining curve, so it is quadratic in the
/// number of spans; it serves curves whose refined knot vector is not in
/// clamped Bezier form.
fn bezier_segments_by_splitting(nurbs: &NurbsCurve) -> Vec<(NurbsCurve, (f64, f64))> {
    let (interior_knots, knots_to_insert) = interior_knots_and_insertions(nurbs);
    if interior_knots.is_empty() {
        return vec![(nurbs.clone(), nurbs.parameter_range())];
    }
    let refined = if knots_to_insert.is_empty() {
        nurbs.clone()
    } else {
        match nurbs.insert_knots(&knots_to_insert) {
            Ok(c) => c,
            Err(_) => nurbs.clone(),
        }
    };

    // Now split at each interior knot.
    let mut segments: Vec<(NurbsCurve, (f64, f64))> = Vec::new();
    let mut remaining = refined;

    for &k in &interior_knots {
        let (r_min, _r_max) = remaining.parameter_range();
        match remaining.split_at(k) {
            Ok((left, right)) => {
                let left_range = left.parameter_range();
                segments.push((left, (r_min.max(left_range.0), k.min(left_range.1))));
                remaining = right;
            }
            Err(_) => {
                // If split fails (e.g. knot at boundary), skip.
            }
        }
    }

    // Add the last segment.
    let rem_range = remaining.parameter_range();
    segments.push((remaining, rem_range));

    segments
}

#[cfg(test)]
#[path = "surface_curve_bezier_tests.rs"]
mod tests;
