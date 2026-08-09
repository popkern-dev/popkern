//! Exact degree/knot compatibility operations for prepared exact-boundary
//! curves: target-degree elevation, union interior knots, knot insertion,
//! budget checks, and prepared-curve verification.

use popkern_core::geometry::NurbsCurve;
use popkern_core::geometry::curve::CurveEval;

use super::{
    ExactBoundaryCurveBudget, ExactBoundaryCurveDeclineReason, ExactBoundaryCurveDiagnostic,
    KNOT_EPS, PreparedSourceCurve, SUPPORT_RIBBON_REFINEMENT_KNOTS, VERIFY_SAMPLES,
    VERIFY_TOLERANCE,
};

/// Domain-relative window at each parameter-domain end inside which an
/// interior knot is a TERMINAL SLIVER, not a feature. Measured
/// (PK-255, cube [0.05,0.1,0.2] at the (1,1,1) corner): a fillet-end
/// source curve carried a full-multiplicity break 2.98e-8 inside its
/// end, the knot union propagated it into every prepared curve, and
/// the derivative-row lever arm (first knot span / degree) collapsed
/// to 3.3e-9 — row 1 fell onto row 0 and every sector's seam G1 died
/// at 0.72–0.99 while the mirrored corner (sliver-free sources)
/// served at e-8. Geometry below this window is sub-tolerance by the
/// kernel's own default (1e-6) and must not shape the solve basis.
const END_SLIVER_WINDOW: f64 = 1.0e-6;

/// Geometric tolerance for removing terminal-sliver knots: above the
/// measured sliver class (~1e-9 3-D extent), below the default kernel
/// tolerance, so a genuine feature declines removal and is kept.
const END_SLIVER_REMOVAL_TOLERANCE: f64 = 1.0e-7;

/// Maximum 3-D displacement accepted when a terminal sliver is
/// trimmed and the retained end control is snapped back onto the
/// original endpoint. A sliver inside [`END_SLIVER_WINDOW`] spans at
/// most ~window × |curve'|, so genuine geometry never reaches this;
/// the pipeline's `verify_prepared_curve` (1e-9) remains the final
/// fidelity gate.
const END_SLIVER_SNAP_TOLERANCE: f64 = 1.0e-6;

/// Remove interior knots that sit within [`END_SLIVER_WINDOW`] of a
/// parameter-domain end. Smooth slivers are removed exactly via
/// tolerance-gated knot removal; kinked sub-tolerance terminal hooks
/// (removal declines — shape-preservation is impossible and unwanted)
/// are TRIMMED at the sliver knot with the retained end control
/// snapped onto the original endpoint, discarding only sub-tolerance
/// geometry. A sliver that survives both stays this curve's own:
/// [`compatible_interior_knots`] refuses to propagate it to sibling
/// curves. The prepared curve is still verified against its source by
/// `verify_prepared_curve`.
pub(super) fn remove_end_sliver_knots(curve: NurbsCurve) -> (NurbsCurve, bool) {
    let mut curve = curve;
    let mut repaired = false;
    loop {
        let (start, end) = curve.parameter_range();
        let window = (end - start).abs() * END_SLIVER_WINDOW;
        let Some((index, multiplicity)) = terminal_sliver_knot(&curve, start, end, window) else {
            return (curve, repaired);
        };
        let value = curve.knots[index];
        let (candidate, removed) =
            curve.remove_knot_times(index, multiplicity, END_SLIVER_REMOVAL_TOLERANCE);
        if removed == multiplicity {
            curve = candidate;
            repaired = true;
            continue;
        }
        let Some(trimmed) = trim_terminal_sliver(&curve, value, start, end) else {
            return if removed > 0 {
                (candidate, true)
            } else {
                (curve, repaired)
            };
        };
        curve = trimmed;
        repaired = true;
    }
}

/// Split the curve at a terminal sliver knot, drop the sub-tolerance
/// terminal piece, snap the retained end control onto the original
/// endpoint, and renormalize the domain to [0, 1]. Returns `None`
/// (keep the original) when the split fails or the snap would move
/// the endpoint by more than [`END_SLIVER_SNAP_TOLERANCE`].
fn trim_terminal_sliver(
    curve: &NurbsCurve,
    sliver_value: f64,
    start: f64,
    end: f64,
) -> Option<NurbsCurve> {
    let end_side = (sliver_value - end).abs() < (sliver_value - start).abs();
    let (left, right) = curve.split_at(sliver_value).ok()?;
    let mut kept = if end_side { left } else { right };
    let original_endpoint = if end_side {
        curve.point_at(end)
    } else {
        curve.point_at(start)
    };
    let snap_target_index = if end_side {
        kept.control_points.len() - 1
    } else {
        0
    };
    let snap_distance = kept.control_points[snap_target_index].distance_to(&original_endpoint);
    if snap_distance > END_SLIVER_SNAP_TOLERANCE {
        return None;
    }
    kept.control_points[snap_target_index] = original_endpoint;
    let (kept_start, kept_end) = kept.parameter_range();
    let span = kept_end - kept_start;
    if span.abs() <= f64::EPSILON {
        return None;
    }
    kept.knots = kept
        .knots
        .iter()
        .map(|knot| (*knot - kept_start) / span)
        .collect();
    Some(kept)
}

fn terminal_sliver_knot(
    curve: &NurbsCurve,
    start: f64,
    end: f64,
    window: f64,
) -> Option<(usize, usize)> {
    let (index, value) = curve.knots.iter().enumerate().find_map(|(index, &knot)| {
        let interior = knot > start + KNOT_EPS && knot < end - KNOT_EPS;
        let near_end = knot <= start + window || knot >= end - window;
        (interior && near_end).then_some((index, knot))
    })?;
    let multiplicity = curve
        .knots
        .iter()
        .filter(|&&knot| (knot - value).abs() <= KNOT_EPS)
        .count();
    Some((index, multiplicity))
}

pub(super) fn elevate_to_target(
    mut curve: PreparedSourceCurve,
    budget: ExactBoundaryCurveBudget,
) -> Result<PreparedSourceCurve, ExactBoundaryCurveDiagnostic> {
    if curve.prepared.degree > budget.target_degree {
        check_budget(&curve, budget)?;
        return Err(ExactBoundaryCurveDiagnostic::new(
            ExactBoundaryCurveDeclineReason::DegreeCompatibilityDeclined,
            &curve,
            None,
        ));
    }
    if curve.prepared.degree < budget.target_degree {
        curve.prepared = curve
            .prepared
            .elevate_degree_to(budget.target_degree)
            .map_err(|_| {
                ExactBoundaryCurveDiagnostic::new(
                    ExactBoundaryCurveDeclineReason::DegreeCompatibilityDeclined,
                    &curve,
                    None,
                )
            })?;
    }
    check_budget(&curve, budget)?;
    Ok(curve)
}

pub(super) fn compatible_interior_knots(curves: &[PreparedSourceCurve]) -> Vec<f64> {
    let mut knots = Vec::new();
    for curve in curves {
        let (start, end) = curve.prepared.parameter_range();
        // Terminal-sliver knots a curve kept (removal declined) stay
        // that curve's own: propagating them would rebuild the nano
        // first/last spans in every sibling basis.
        let window = (end - start).abs() * END_SLIVER_WINDOW;
        for &knot in &curve.prepared.knots {
            if knot > start + window
                && knot < end - window
                && !knots
                    .iter()
                    .any(|existing| f64::abs(*existing - knot) <= KNOT_EPS)
            {
                knots.push(knot);
            }
        }
    }
    knots.extend(SUPPORT_RIBBON_REFINEMENT_KNOTS);
    knots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    knots
}

pub(super) fn insert_missing_knots(
    curve: &mut PreparedSourceCurve,
    union_knots: &[f64],
    budget: ExactBoundaryCurveBudget,
) -> Result<(), ExactBoundaryCurveDiagnostic> {
    for &knot in union_knots {
        let required = union_knots
            .iter()
            .filter(|required| f64::abs(**required - knot) <= KNOT_EPS)
            .count();
        let existing = curve
            .prepared
            .knots
            .iter()
            .filter(|existing| f64::abs(**existing - knot) <= KNOT_EPS)
            .count();
        if existing < required {
            curve.prepared = curve.prepared.insert_knot(knot).map_err(|_| {
                ExactBoundaryCurveDiagnostic::new(
                    ExactBoundaryCurveDeclineReason::KnotCompatibilityDeclined,
                    curve,
                    None,
                )
            })?;
            check_budget(curve, budget)?;
        }
    }
    Ok(())
}

pub(super) fn check_budget(
    curve: &PreparedSourceCurve,
    budget: ExactBoundaryCurveBudget,
) -> Result<(), ExactBoundaryCurveDiagnostic> {
    if curve.prepared.control_points.len() > budget.max_control_points
        || curve.prepared.knots.len() > budget.max_knots
    {
        Err(ExactBoundaryCurveDiagnostic::new(
            ExactBoundaryCurveDeclineReason::KnotBudgetExceeded,
            curve,
            None,
        ))
    } else {
        Ok(())
    }
}

/// Fidelity gate for curves that underwent terminal-sliver repair:
/// the repair intentionally discards sub-tolerance geometry (snap ≤
/// [`END_SLIVER_SNAP_TOLERANCE`], reparameterization drift ≤ window ×
/// curve length), so exactness verification at `VERIFY_TOLERANCE` is
/// impossible by design. This bound is the kernel's default geometric
/// tolerance — anything worse is a genuine preparation failure.
const SLIVER_REPAIR_VERIFY_TOLERANCE: f64 = 1.0e-6;

pub(super) fn verify_prepared_curve(
    curve: &PreparedSourceCurve,
) -> Result<(), ExactBoundaryCurveDiagnostic> {
    let (start, end) = curve.source_range;
    let mut worst = 0.0;
    for i in 0..VERIFY_SAMPLES {
        let normalized = i as f64 / (VERIFY_SAMPLES - 1) as f64;
        let t = start + (end - start) * normalized;
        let residual = curve
            .source
            .point_at(t)
            .distance_to(&curve.prepared.point_at(normalized));
        if residual > worst {
            worst = residual;
        }
    }
    let tolerance = if curve.sliver_repaired {
        SLIVER_REPAIR_VERIFY_TOLERANCE
    } else {
        VERIFY_TOLERANCE
    };
    if worst > tolerance {
        Err(ExactBoundaryCurveDiagnostic::new(
            ExactBoundaryCurveDeclineReason::BoundaryVerificationFailed,
            curve,
            Some(worst),
        ))
    } else {
        Ok(())
    }
}
