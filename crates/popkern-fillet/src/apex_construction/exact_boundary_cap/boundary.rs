//! Exact-boundary curve preparation for apex cap sectors.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::healing::deviation::edge_surface_deviation;
use popkern_core::topology::{BRep, EdgeHandle};

use super::sector::{ExactBoundaryCapSectorInput, ExactBoundaryOuterSeam, ExactBoundaryOuterSide};

const VERIFY_SAMPLES: usize = 9;
const VERIFY_TOLERANCE: f64 = 1.0e-9;
const KNOT_EPS: f64 = 1.0e-10;
const TOPOLOGY_ENDPOINT_TOLERANCE: f64 = 1.0e-7;
const SUPPORT_SOURCE_SAMPLES: usize = 17;
const SUPPORT_SOURCE_TOLERANCE: f64 = 1.0e-6;

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
                sector_name: self.sector_name,
                edge: Some(self.edge),
                support_half_edge: None,
                support_face: None,
                residual: self.residual,
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
    if let Some(support_curve) = validated_support_pcurve_source_curve(brep, seam, &raw_curve) {
        return prepared_source_from_curve(brep, sector, seam, support_curve);
    }
    prepared_source_from_curve(brep, sector, seam, raw_curve)
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

fn validated_support_pcurve_source_curve(
    brep: &BRep,
    seam: &ExactBoundaryOuterSeam,
    edge_curve: &NurbsCurve,
) -> Option<NurbsCurve> {
    let support_curve = support_pcurve_source_curve(brep, seam)?;
    let (origin, destination) = cap_half_edge_vertex_span(brep, seam.cap_half_edge)?;
    if endpoint_orientation_gap(&support_curve, origin, destination) > TOPOLOGY_ENDPOINT_TOLERANCE
        && endpoint_orientation_gap(&reverse_nurbs_curve(&support_curve), origin, destination)
            > TOPOLOGY_ENDPOINT_TOLERANCE
    {
        return None;
    }
    let support_surface = support_surface_for_seam(brep, seam)?;
    let support_range = support_curve.parameter_range();
    let support_deviation = edge_surface_deviation(
        &support_curve,
        support_range,
        support_surface,
        SUPPORT_SOURCE_SAMPLES,
    )
    .max_deviation;
    if support_deviation > SUPPORT_SOURCE_TOLERANCE {
        return None;
    }
    let edge_range = edge_curve.parameter_range();
    let edge_deviation = edge_surface_deviation(
        edge_curve,
        edge_range,
        support_surface,
        SUPPORT_SOURCE_SAMPLES,
    )
    .max_deviation;
    let cap_pcurve_missing = brep
        .half_edges
        .get(seam.cap_half_edge)
        .is_some_and(|half_edge| half_edge.pcurve.is_none());
    (edge_deviation > support_deviation.max(SUPPORT_SOURCE_TOLERANCE) || cap_pcurve_missing)
        .then_some(support_curve)
}

fn support_pcurve_source_curve(brep: &BRep, seam: &ExactBoundaryOuterSeam) -> Option<NurbsCurve> {
    let support_half_edge = brep.half_edges.get(seam.support_half_edge)?;
    let pcurve = support_half_edge.pcurve.as_ref()?;
    let support_surface = support_surface_for_seam(brep, seam)?;
    let (t0, t1) = pcurve.parameter_range();
    let mut points = Vec::with_capacity(SUPPORT_SOURCE_SAMPLES);
    let mut params = Vec::with_capacity(SUPPORT_SOURCE_SAMPLES);
    for index in 0..SUPPORT_SOURCE_SAMPLES {
        let parameter = index as f64 / (SUPPORT_SOURCE_SAMPLES - 1) as f64;
        let pcurve_t = t0 + (t1 - t0) * parameter;
        let uv = pcurve.point_at(pcurve_t);
        points.push(support_surface.point_at(uv.x, uv.y));
        params.push(parameter);
    }
    let degree = 9.min(points.len().saturating_sub(1));
    popkern_core::geometry::fitting::fit_curve(&points, &params, degree, points.len()).ok()
}

fn support_surface_for_seam<'a>(
    brep: &'a BRep,
    seam: &ExactBoundaryOuterSeam,
) -> Option<&'a dyn SurfaceEval> {
    let support_half_edge = brep.half_edges.get(seam.support_half_edge)?;
    let support_face = brep.loops.get(support_half_edge.loop_)?.face;
    Some(&brep.faces.get(support_face)?.surface.as_ref()?.geometry)
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

fn cap_half_edge_vertex_span(
    brep: &BRep,
    half_edge: popkern_core::topology::HalfEdgeHandle,
) -> Option<(Point3, Point3)> {
    let half_edge_data = brep.half_edges.get(half_edge)?;
    let twin = brep.half_edges.get(half_edge_data.twin)?;
    let origin = brep.vertices.get(half_edge_data.origin)?.position;
    let destination = brep.vertices.get(twin.origin)?.position;
    Some((origin, destination))
}

fn endpoint_orientation_gap(curve: &NurbsCurve, origin: Point3, destination: Point3) -> f64 {
    let (start, end) = curve.parameter_range();
    curve
        .point_at(start)
        .distance_to(&origin)
        .max(curve.point_at(end).distance_to(&destination))
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

fn elevate_to_target(
    mut curve: PreparedSourceCurve,
    budget: ExactBoundaryCurveBudget,
) -> Result<PreparedSourceCurve, ExactBoundaryCurveDiagnostic> {
    if curve.prepared.degree > budget.target_degree {
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

fn compatible_interior_knots(curves: &[PreparedSourceCurve]) -> Vec<f64> {
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
    knots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    knots
}

fn insert_missing_knots(
    curve: &mut PreparedSourceCurve,
    union_knots: &[f64],
    budget: ExactBoundaryCurveBudget,
) -> Result<(), ExactBoundaryCurveDiagnostic> {
    for &knot in union_knots {
        if !curve
            .prepared
            .knots
            .iter()
            .any(|existing| f64::abs(*existing - knot) <= KNOT_EPS)
        {
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

fn check_budget(
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

fn verify_prepared_curve(curve: &PreparedSourceCurve) -> Result<(), ExactBoundaryCurveDiagnostic> {
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

fn reverse_nurbs_curve(curve: &NurbsCurve) -> NurbsCurve {
    let first = *curve.knots.first().unwrap_or(&0.0);
    let last = *curve.knots.last().unwrap_or(&1.0);
    let mut knots: Vec<f64> = curve
        .knots
        .iter()
        .rev()
        .map(|k| first + last - *k)
        .collect();
    knots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mut control_points: Vec<Point3> = curve.control_points.clone();
    control_points.reverse();
    let mut weights = curve.weights.clone();
    weights.reverse();
    NurbsCurve {
        degree: curve.degree,
        knots,
        control_points,
        weights,
    }
}
