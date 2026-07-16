//! Exact degree/knot compatibility operations for prepared exact-boundary
//! curves: target-degree elevation, union interior knots, knot insertion,
//! budget checks, and prepared-curve verification.

use popkern_core::geometry::curve::CurveEval;

use super::{
    ExactBoundaryCurveBudget, ExactBoundaryCurveDeclineReason, ExactBoundaryCurveDiagnostic,
    KNOT_EPS, PreparedSourceCurve, SUPPORT_RIBBON_REFINEMENT_KNOTS, VERIFY_SAMPLES,
    VERIFY_TOLERANCE,
};

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
        for &knot in &curve.prepared.knots {
            if knot > start + KNOT_EPS
                && knot < end - KNOT_EPS
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
    if worst > VERIFY_TOLERANCE {
        Err(ExactBoundaryCurveDiagnostic::new(
            ExactBoundaryCurveDeclineReason::BoundaryVerificationFailed,
            curve,
            Some(worst),
        ))
    } else {
        Ok(())
    }
}
