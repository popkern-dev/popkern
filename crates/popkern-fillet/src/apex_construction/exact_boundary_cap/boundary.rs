use popkern_core::geometry::NurbsCurve;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::topology::{BRep, EdgeHandle};

#[path = "boundary_compatibility.rs"]
mod boundary_compatibility;
#[path = "boundary_curve_orientation.rs"]
mod boundary_curve_orientation;
#[path = "boundary_support_source.rs"]
mod boundary_support_source;

use super::sector::{ExactBoundaryCapSectorInput, ExactBoundaryOuterSeam, ExactBoundaryOuterSide};
use boundary_compatibility::{
    check_budget, compatible_interior_knots, elevate_to_target, insert_missing_knots,
    verify_prepared_curve,
};
use boundary_curve_orientation::{
    cap_half_edge_vertex_span, endpoint_orientation_gap, reverse_nurbs_curve,
};
use boundary_support_source::{
    raw_curve_needs_support_authority, validated_projected_support_source_curve,
    validated_support_pcurve_source_curve,
};

const VERIFY_SAMPLES: usize = 9;
const VERIFY_TOLERANCE: f64 = 1.0e-9;
const KNOT_EPS: f64 = 1.0e-10;
const TOPOLOGY_ENDPOINT_TOLERANCE: f64 = 1.0e-7;
#[rustfmt::skip]
const SUPPORT_RIBBON_REFINEMENT_KNOTS: [f64; 28] = [0.125, 0.125, 0.125, 0.125, 0.250, 0.250, 0.250, 0.250, 0.375, 0.375, 0.375, 0.375, 0.500, 0.500, 0.500, 0.500, 0.625, 0.625, 0.625, 0.625, 0.750, 0.750, 0.750, 0.750, 0.875, 0.875, 0.875, 0.875];

/// Budget for exact boundary compatibility operations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct ExactBoundaryCurveBudget {
    /// Required common degree. Curves already above this degree decline.
    pub target_degree: usize,
    /// Maximum allowed control points after exact operations.
    pub max_control_points: usize,
    /// Maximum allowed knots after exact operations.
    pub max_knots: usize,
}

impl Default for ExactBoundaryCurveBudget {
    fn default() -> Self {
        Self {
            target_degree: 9,
            max_control_points: 128,
            max_knots: 140,
        }
    }
}

/// Prepared exact boundary curve for one immutable outer seam.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct PreparedExactBoundaryCurve {
    /// Planner-order sector index.
    pub sector_index: usize,
    /// Stable sector name.
    pub sector_name: String,
    /// Sector-local side.
    pub side: ExactBoundaryOuterSide,
    /// Source BRep edge.
    pub edge: EdgeHandle,
    /// Oriented and compatibility-prepared curve.
    pub curve: NurbsCurve,
}

/// Prepared exact boundary curves for all extracted sectors.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct PreparedExactBoundaryCurves {
    /// One prepared curve per immutable outer seam.
    pub curves: Vec<PreparedExactBoundaryCurve>,
}

/// Stable decline reason for exact boundary curve preparation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ExactBoundaryCurveDeclineReason {
    /// Source edge curve was missing or stale.
    BoundaryMissingCurve,
    /// Exact degree elevation could not reach the target degree.
    DegreeCompatibilityDeclined,
    /// Knot insertion failed.
    KnotCompatibilityDeclined,
    /// Exact operation exceeded configured budget.
    KnotBudgetExceeded,
    /// Prepared curve no longer evaluates to the source curve.
    BoundaryVerificationFailed,
}

impl ExactBoundaryCurveDeclineReason {
    /// Stable diagnostic code for tests and downstream APIs.
    pub(crate) fn code(self) -> &'static str {
        match self {
            Self::BoundaryMissingCurve => "exact_boundary_cap.boundary_missing_curve",
            Self::DegreeCompatibilityDeclined => "exact_boundary_cap.knot_compatibility_declined",
            Self::KnotCompatibilityDeclined => "exact_boundary_cap.knot_compatibility_declined",
            Self::KnotBudgetExceeded => "exact_boundary_cap.knot_budget_exceeded",
            Self::BoundaryVerificationFailed => "exact_boundary_cap.boundary_orientation_conflict",
        }
    }
}

/// Structured exact-boundary curve preparation diagnostic.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ExactBoundaryCurveDiagnostic {
    /// Stable decline reason.
    pub reason: ExactBoundaryCurveDeclineReason,
    /// Stable diagnostic code.
    pub code: &'static str,
    /// Sector index being prepared.
    pub sector_index: usize,
    /// Stable sector name.
    pub sector_name: String,
    /// Sector-local side.
    pub side: ExactBoundaryOuterSide,
    /// Source BRep edge.
    pub edge: EdgeHandle,
    /// Measured verification residual, when applicable.
    pub residual: Option<f64>,
}

impl ExactBoundaryCurveDiagnostic {
    fn new(
        reason: ExactBoundaryCurveDeclineReason,
        source: &PreparedSourceCurve,
        residual: Option<f64>,
    ) -> Self {
        Self {
            reason,
            code: reason.code(),
            sector_index: source.sector_index,
            sector_name: source.sector_name.clone(),
            side: source.side,
            edge: source.edge,
            residual,
        }
    }

    pub(crate) fn into_apex_diagnostic(self) -> crate::types::ApexConstructionDiagnostic {
        crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(
            crate::types::ExactBoundaryCapDiagnostic {
                stage: "exact_boundary_cap.boundary".to_string(),
                code: self.code.to_string(),
                detail_code: None,
                sector_name: self.sector_name,
                edge: Some(self.edge),
                support_half_edge: None,
                support_face: None,
                residual: self.residual,
                outer_side: Some(format!("{:?}", self.side)),
                sample_parameter: None,
                tolerance: None,
                shared_fact_ids: Vec::new(),
            },
        )
    }
}

struct PreparedSourceCurve {
    sector_index: usize,
    sector_name: String,
    side: ExactBoundaryOuterSide,
    edge: EdgeHandle,
    source: NurbsCurve,
    source_range: (f64, f64),
    prepared: NurbsCurve,
}

/// Prepare exact boundary curves using only orientation normalization, degree
/// elevation, and knot insertion.
pub(crate) fn prepare_exact_boundary_curves(
    brep: &BRep,
    sectors: &[ExactBoundaryCapSectorInput],
    budget: ExactBoundaryCurveBudget,
) -> Result<PreparedExactBoundaryCurves, ExactBoundaryCurveDiagnostic> {
    let mut curves = Vec::with_capacity(sectors.len() * 2);
    for sector in sectors {
        for seam in &sector.outer_seams {
            let source = oriented_source_curve(brep, sector, seam)?;
            let prepared = elevate_to_target(source, budget)?;
            curves.push(prepared);
        }
    }

    let union_knots = compatible_interior_knots(&curves);
    for curve in &mut curves {
        insert_missing_knots(curve, &union_knots, budget)?;
        check_budget(curve, budget)?;
        verify_prepared_curve(curve)?;
    }

    Ok(PreparedExactBoundaryCurves {
        curves: curves
            .into_iter()
            .map(|source| PreparedExactBoundaryCurve {
                sector_index: source.sector_index,
                sector_name: source.sector_name,
                side: source.side,
                edge: source.edge,
                curve: source.prepared,
            })
            .collect(),
    })
}

fn oriented_source_curve(
    brep: &BRep,
    sector: &ExactBoundaryCapSectorInput,
    seam: &ExactBoundaryOuterSeam,
) -> Result<PreparedSourceCurve, ExactBoundaryCurveDiagnostic> {
    let side = seam.side;
    let edge_data = brep
        .edges
        .get(seam.edge)
        .ok_or_else(|| ExactBoundaryCurveDiagnostic {
            reason: ExactBoundaryCurveDeclineReason::BoundaryMissingCurve,
            code: ExactBoundaryCurveDeclineReason::BoundaryMissingCurve.code(),
            sector_index: sector.sector_index,
            sector_name: sector.sector_name.clone(),
            side,
            edge: seam.edge,
            residual: None,
        })?;
    let Some(curve_on_edge) = edge_data.curve.as_ref() else {
        return Err(ExactBoundaryCurveDiagnostic {
            reason: ExactBoundaryCurveDeclineReason::BoundaryMissingCurve,
            code: ExactBoundaryCurveDeclineReason::BoundaryMissingCurve.code(),
            sector_index: sector.sector_index,
            sector_name: sector.sector_name.clone(),
            side,
            edge: seam.edge,
            residual: None,
        });
    };
    let raw_curve = curve_on_edge
        .geometry
        .to_nurbs_in_range((curve_on_edge.t_start, curve_on_edge.t_end));
    // The support pcurve is the boundary authority whenever the stored
    // edge curve is measurably off its support surface — including exact
    // rational arcs, whose termination construction can land ~1e-5 off
    // the fillet cylinder and would otherwise bake that error into the
    // cap boundary row (PK-124 marginal FilletEnd seams).
    let needs_support_authority = raw_curve_needs_support_authority(brep, seam, &raw_curve);
    if (!has_non_unit_weights(&raw_curve) || needs_support_authority)
        && let Some(support_curve) = validated_support_pcurve_source_curve(brep, seam, &raw_curve)
        && let Ok(prepared) = prepared_source_from_curve(brep, sector, seam, support_curve)
    {
        if needs_support_authority && crate::direct_apex_trace::enabled() {
            eprintln!(
                "[direct-apex-trace exact-boundary/boundary] seam edge {:?}: off-support raw \
                 curve replaced by support-pcurve authority",
                seam.edge
            );
        }
        return Ok(prepared);
    }
    match prepared_source_from_curve(brep, sector, seam, raw_curve) {
        Ok(prepared) => Ok(prepared),
        Err(diagnostic)
            if diagnostic.reason == ExactBoundaryCurveDeclineReason::BoundaryVerificationFailed =>
        {
            if let Some(topology_curve) = topology_span_source_curve(brep, seam) {
                if let Some(support_curve) =
                    validated_projected_support_source_curve(brep, seam, &topology_curve)
                    && let Ok(prepared) =
                        prepared_source_from_curve(brep, sector, seam, support_curve)
                {
                    return Ok(prepared);
                }
                return prepared_source_from_curve(brep, sector, seam, topology_curve);
            }
            Err(diagnostic)
        }
        Err(diagnostic) => Err(diagnostic),
    }
}

fn has_non_unit_weights(curve: &NurbsCurve) -> bool {
    curve.weights.len() == curve.control_points.len()
        && curve
            .weights
            .iter()
            .any(|weight| (*weight - 1.0).abs() > 1.0e-9)
}

fn topology_span_source_curve(brep: &BRep, seam: &ExactBoundaryOuterSeam) -> Option<NurbsCurve> {
    let (origin, destination) = cap_half_edge_vertex_span(brep, seam.cap_half_edge)?;
    if origin.distance_to(&destination) <= TOPOLOGY_ENDPOINT_TOLERANCE {
        return None;
    }
    Some(NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![origin, destination],
        weights: vec![1.0, 1.0],
    })
}

fn prepared_source_from_curve(
    brep: &BRep,
    sector: &ExactBoundaryCapSectorInput,
    seam: &ExactBoundaryOuterSeam,
    raw_curve: NurbsCurve,
) -> Result<PreparedSourceCurve, ExactBoundaryCurveDiagnostic> {
    let side = seam.side;
    let curve = orient_curve_to_cap_half_edge(brep, sector, seam, raw_curve)?;
    let source_range = curve.parameter_range();
    let prepared = normalize_curve_parameter_domain(&curve).ok_or_else(|| {
        ExactBoundaryCurveDiagnostic::new(
            ExactBoundaryCurveDeclineReason::BoundaryVerificationFailed,
            &PreparedSourceCurve {
                sector_index: sector.sector_index,
                sector_name: sector.sector_name.clone(),
                side,
                edge: seam.edge,
                source: curve.clone(),
                source_range,
                prepared: curve.clone(),
            },
            None,
        )
    })?;
    Ok(PreparedSourceCurve {
        sector_index: sector.sector_index,
        sector_name: sector.sector_name.clone(),
        side,
        edge: seam.edge,
        source: curve,
        source_range,
        prepared,
    })
}

fn orient_curve_to_cap_half_edge(
    brep: &BRep,
    sector: &ExactBoundaryCapSectorInput,
    seam: &ExactBoundaryOuterSeam,
    curve: NurbsCurve,
) -> Result<NurbsCurve, ExactBoundaryCurveDiagnostic> {
    let Some((origin, destination)) = cap_half_edge_vertex_span(brep, seam.cap_half_edge) else {
        return Err(ExactBoundaryCurveDiagnostic::new(
            ExactBoundaryCurveDeclineReason::BoundaryVerificationFailed,
            &diagnostic_source_curve(sector, seam, curve, None),
            None,
        ));
    };
    let direct_gap = endpoint_orientation_gap(&curve, origin, destination);
    if direct_gap <= TOPOLOGY_ENDPOINT_TOLERANCE {
        return Ok(curve);
    }

    let reversed = reverse_nurbs_curve(&curve);
    let reverse_gap = endpoint_orientation_gap(&reversed, origin, destination);
    if reverse_gap <= TOPOLOGY_ENDPOINT_TOLERANCE || reverse_gap < direct_gap {
        return Ok(reversed);
    }

    Err(ExactBoundaryCurveDiagnostic::new(
        ExactBoundaryCurveDeclineReason::BoundaryVerificationFailed,
        &diagnostic_source_curve(sector, seam, curve, None),
        Some(direct_gap.min(reverse_gap)),
    ))
}

fn diagnostic_source_curve(
    sector: &ExactBoundaryCapSectorInput,
    seam: &ExactBoundaryOuterSeam,
    curve: NurbsCurve,
    source_range: Option<(f64, f64)>,
) -> PreparedSourceCurve {
    let source_range = source_range.unwrap_or_else(|| curve.parameter_range());
    PreparedSourceCurve {
        sector_index: sector.sector_index,
        sector_name: sector.sector_name.clone(),
        side: seam.side,
        edge: seam.edge,
        source: curve.clone(),
        source_range,
        prepared: curve,
    }
}

fn normalize_curve_parameter_domain(curve: &NurbsCurve) -> Option<NurbsCurve> {
    let (start, end) = curve.parameter_range();
    let span = end - start;
    if span.abs() <= f64::EPSILON {
        return None;
    }
    Some(NurbsCurve {
        degree: curve.degree,
        knots: curve
            .knots
            .iter()
            .map(|knot| (*knot - start) / span)
            .collect(),
        control_points: curve.control_points.clone(),
        weights: curve.weights.clone(),
    })
}
