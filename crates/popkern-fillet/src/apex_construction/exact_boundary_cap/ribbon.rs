//! Support derivative ribbon preparation for exact-boundary apex caps.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::healing::deviation::edge_surface_deviation;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

use super::boundary::{PreparedExactBoundaryCurve, PreparedExactBoundaryCurves};
#[cfg(test)]
pub(super) use super::ribbon_projection::compute_support_ribbon_sample;
use super::ribbon_projection::{
    choose_support_pcurve_orientation, compute_support_ribbon_samples, local_sector_width,
    projected_support_pcurve,
};
use super::sector::{ExactBoundaryCapSectorInput, ExactBoundaryOuterSide};
use super::support::{ExactBoundarySupportSurfaceFamily, support_surface_family};

#[path = "ribbon_surface_identity.rs"]
mod ribbon_surface_identity;
use ribbon_surface_identity::geometrically_unique_surface;

const DEFAULT_SAMPLE_COUNT: usize = 5;

/// Budget and local tolerances for support derivative ribbon preparation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct ExactBoundaryRibbonBudget {
    /// Number of samples per immutable support seam.
    pub sample_count: usize,
    /// Maximum accepted distance between edge curve and support pcurve image.
    pub max_support_residual: f64,
    /// Minimum support normal length before normalization.
    pub min_normal_length: f64,
    /// Minimum along-boundary tangent length after tangent-plane projection.
    pub min_tangent_length: f64,
    /// Minimum cross-boundary direction length.
    pub min_cross_length: f64,
    /// Maximum derivative vector length used by this preparatory slice.
    pub max_derivative_scale: f64,
}

impl Default for ExactBoundaryRibbonBudget {
    fn default() -> Self {
        Self {
            sample_count: DEFAULT_SAMPLE_COUNT,
            max_support_residual: 1.0e-6,
            min_normal_length: 1.0e-12,
            min_tangent_length: 1.0e-12,
            min_cross_length: 1.0e-12,
            max_derivative_scale: 1.0,
        }
    }
}

/// Prepared support derivative ribbons for all exact-boundary outer seams.
#[derive(Debug, Clone)]
pub(crate) struct PreparedExactBoundaryRibbons {
    /// One ribbon per immutable outer seam.
    pub ribbons: Vec<PreparedExactBoundaryRibbon>,
}

/// Support derivative ribbon for one immutable outer seam.
#[derive(Debug, Clone)]
pub(crate) struct PreparedExactBoundaryRibbon {
    /// Planner-order sector index.
    pub sector_index: usize,
    /// Stable sector name.
    pub sector_name: String,
    /// Sector-local side.
    pub side: ExactBoundaryOuterSide,
    /// Source BRep edge.
    pub edge: EdgeHandle,
    /// Support-side half-edge that supplied the pcurve evidence.
    pub support_half_edge: HalfEdgeHandle,
    /// Support face sampled by the pcurve.
    pub support_face: FaceHandle,
    /// Support-surface family used for support-ribbon evaluation.
    pub support_family: ExactBoundarySupportSurfaceFamily,
    /// Exact support surface geometry sampled for this seam.
    pub support_surface: SurfaceGeometry,
    /// True when support pcurve orientation is opposite the prepared cap curve.
    pub support_parameter_reversed: bool,
    /// True when support half-edge orientation matches the edge curve parameter direction.
    pub support_half_edge_curve_forward: bool,
    /// Worst sampled 3-D residual between exact edge curve and pcurve image.
    pub max_support_residual: f64,
    /// Sampled derivative data along the seam.
    pub samples: Vec<ExactBoundaryRibbonSample>,
}

/// One sampled support Hermite datum along an exact boundary seam.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct ExactBoundaryRibbonSample {
    /// Normalized prepared-boundary parameter.
    pub parameter: f64,
    /// Support pcurve UV sample.
    pub support_uv: (f64, f64),
    /// Exact boundary point from the prepared edge curve.
    pub boundary_point: Point3,
    /// Support surface point from pcurve evidence.
    pub support_point: Point3,
    /// Unit along-boundary tangent projected into the support tangent plane.
    pub boundary_tangent: Vec3,
    /// Unit support normal.
    pub support_normal: Vec3,
    /// Cross-boundary derivative vector in the support tangent plane.
    pub cross_derivative: Vec3,
}

/// Stable decline reason for support derivative ribbon preparation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ExactBoundaryRibbonDeclineReason {
    /// Prepared boundary curve for the sector seam was missing.
    BoundaryCurveMissing,
    /// Support face or support surface was missing/stale.
    SupportSurfaceMissing,
    /// Support half-edge pcurve was missing.
    SupportPcurveMissing,
    /// Support pcurve image did not match the exact edge curve.
    SupportPcurveStale,
    /// Surface normal was degenerate.
    DegenerateSupportNormal,
    /// Edge tangent could not define a support tangent-plane direction.
    TangentAmbiguity,
    /// Cross-boundary derivative direction collapsed.
    CrossDerivativeDegenerate,
    /// Support surface family is not accepted for exact-boundary ribbons yet.
    UnsupportedSupportSurface,
}

impl ExactBoundaryRibbonDeclineReason {
    /// Stable diagnostic code for tests, API, and future UI highlighting.
    pub(crate) fn code(self) -> &'static str {
        match self {
            Self::BoundaryCurveMissing => "exact_boundary_cap.ribbon.boundary_curve_missing",
            Self::SupportSurfaceMissing => "exact_boundary_cap.ribbon.support_surface_missing",
            Self::SupportPcurveMissing => "exact_boundary_cap.ribbon.support_pcurve_missing",
            Self::SupportPcurveStale => "exact_boundary_cap.ribbon.support_pcurve_stale",
            Self::DegenerateSupportNormal => "exact_boundary_cap.ribbon.degenerate_support_normal",
            Self::TangentAmbiguity => "exact_boundary_cap.ribbon.tangent_ambiguity",
            Self::CrossDerivativeDegenerate => {
                "exact_boundary_cap.ribbon.cross_derivative_degenerate"
            }
            Self::UnsupportedSupportSurface => {
                "exact_boundary_cap.ribbon.unsupported_support_surface"
            }
        }
    }
}

/// Structured support-ribbon decline.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ExactBoundaryRibbonDiagnostic {
    /// Stable decline reason.
    pub reason: ExactBoundaryRibbonDeclineReason,
    /// Stable diagnostic code.
    pub code: &'static str,
    /// Planner-order sector index.
    pub sector_index: usize,
    /// Stable sector name.
    pub sector_name: String,
    /// Sector-local side.
    pub side: ExactBoundaryOuterSide,
    /// Source BRep edge.
    pub edge: EdgeHandle,
    /// Support half-edge being sampled.
    pub support_half_edge: HalfEdgeHandle,
    /// Support face being sampled.
    pub support_face: FaceHandle,
    /// Measured residual when available.
    pub residual: Option<f64>,
}

/// Input for one support-ribbon sample.
pub(super) struct SupportRibbonSampleRequest<'a> {
    pub boundary: &'a PreparedExactBoundaryCurve,
    pub other_boundary: &'a PreparedExactBoundaryCurve,
    pub surface: &'a SurfaceGeometry,
    pub pcurve: &'a Pcurve,
    pub parameter: f64,
    pub pcurve_reversed: bool,
    pub derivative_scale: f64,
    pub budget: ExactBoundaryRibbonBudget,
}

impl ExactBoundaryRibbonDiagnostic {
    fn new(
        reason: ExactBoundaryRibbonDeclineReason,
        context: &RibbonContext<'_>,
        residual: Option<f64>,
    ) -> Self {
        Self {
            reason,
            code: reason.code(),
            sector_index: context.sector.sector_index,
            sector_name: context.sector.sector_name.clone(),
            side: context.side,
            edge: context.edge,
            support_half_edge: context.support_half_edge,
            support_face: context.support_face,
            residual,
        }
    }

    pub(crate) fn into_apex_diagnostic(self) -> crate::types::ApexConstructionDiagnostic {
        crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(
            crate::types::ExactBoundaryCapDiagnostic {
                stage: "exact_boundary_cap.support_ribbon".to_string(),
                code: self.code.to_string(),
                sector_name: self.sector_name,
                edge: Some(self.edge),
                support_half_edge: Some(self.support_half_edge),
                support_face: Some(self.support_face),
                residual: self.residual,
                shared_fact_ids: Vec::new(),
            },
        )
    }
}

struct RibbonContext<'a> {
    sector: &'a ExactBoundaryCapSectorInput,
    side: ExactBoundaryOuterSide,
    edge: EdgeHandle,
    support_half_edge: HalfEdgeHandle,
    support_face: FaceHandle,
    support_half_edge_curve_forward: bool,
}

/// Prepare support derivative ribbons from exact boundaries and support pcurves.
pub(crate) fn prepare_exact_boundary_ribbons(
    brep: &BRep,
    sectors: &[ExactBoundaryCapSectorInput],
    boundaries: &PreparedExactBoundaryCurves,
    budget: ExactBoundaryRibbonBudget,
) -> Result<PreparedExactBoundaryRibbons, ExactBoundaryRibbonDiagnostic> {
    let mut ribbons = Vec::with_capacity(sectors.len() * 2);
    for sector in sectors {
        for seam in &sector.outer_seams {
            let context = RibbonContext {
                sector,
                side: seam.side,
                edge: seam.edge,
                support_half_edge: seam.support_half_edge,
                support_face: seam.support_face,
                support_half_edge_curve_forward: support_half_edge_curve_forward(
                    brep,
                    seam.support_half_edge,
                )
                .unwrap_or(true),
            };
            let boundary =
                boundary_for(boundaries, sector.sector_index, seam.side).ok_or_else(|| {
                    ExactBoundaryRibbonDiagnostic::new(
                        ExactBoundaryRibbonDeclineReason::BoundaryCurveMissing,
                        &context,
                        None,
                    )
                })?;
            let other_boundary = other_boundary_for(boundaries, sector.sector_index, seam.side)
                .ok_or_else(|| {
                    ExactBoundaryRibbonDiagnostic::new(
                        ExactBoundaryRibbonDeclineReason::BoundaryCurveMissing,
                        &context,
                        None,
                    )
                })?;
            let surface = support_surface(brep, &context).ok_or_else(|| {
                ExactBoundaryRibbonDiagnostic::new(
                    ExactBoundaryRibbonDeclineReason::SupportSurfaceMissing,
                    &context,
                    None,
                )
            })?;
            let pcurve = brep
                .half_edges
                .get(seam.support_half_edge)
                .and_then(|half_edge| half_edge.pcurve.as_ref());

            ribbons.push(prepare_ribbon_for_seam(
                &context,
                boundary,
                other_boundary,
                &surface,
                pcurve,
                budget,
            )?);
        }
    }
    Ok(PreparedExactBoundaryRibbons { ribbons })
}

fn prepare_ribbon_for_seam(
    context: &RibbonContext<'_>,
    boundary: &PreparedExactBoundaryCurve,
    other_boundary: &PreparedExactBoundaryCurve,
    surface: &SurfaceGeometry,
    pcurve: Option<&Pcurve>,
    budget: ExactBoundaryRibbonBudget,
) -> Result<PreparedExactBoundaryRibbon, ExactBoundaryRibbonDiagnostic> {
    let support_family = support_surface_family(surface);
    if matches!(
        support_family,
        ExactBoundarySupportSurfaceFamily::Unsupported
    ) {
        return Err(ExactBoundaryRibbonDiagnostic::new(
            ExactBoundaryRibbonDeclineReason::UnsupportedSupportSurface,
            context,
            None,
        ));
    }
    let projected_pcurve;
    let (pcurve, orientation) = match pcurve {
        Some(existing) => {
            let orientation =
                choose_support_pcurve_orientation(boundary, surface, existing, budget);
            if orientation.max_residual <= budget.max_support_residual {
                (existing, orientation)
            } else {
                projected_pcurve =
                    projected_support_pcurve(boundary, surface, budget).ok_or_else(|| {
                        ExactBoundaryRibbonDiagnostic::new(
                            ExactBoundaryRibbonDeclineReason::SupportPcurveStale,
                            context,
                            Some(orientation.max_residual),
                        )
                    })?;
                let projected_orientation =
                    choose_support_pcurve_orientation(boundary, surface, &projected_pcurve, budget);
                if projected_orientation.max_residual > budget.max_support_residual {
                    return Err(ExactBoundaryRibbonDiagnostic::new(
                        ExactBoundaryRibbonDeclineReason::SupportPcurveStale,
                        context,
                        Some(projected_orientation.max_residual),
                    ));
                }
                (&projected_pcurve, projected_orientation)
            }
        }
        None => {
            projected_pcurve =
                projected_support_pcurve(boundary, surface, budget).ok_or_else(|| {
                    ExactBoundaryRibbonDiagnostic::new(
                        ExactBoundaryRibbonDeclineReason::SupportPcurveMissing,
                        context,
                        None,
                    )
                })?;
            let projected_orientation =
                choose_support_pcurve_orientation(boundary, surface, &projected_pcurve, budget);
            if projected_orientation.max_residual > budget.max_support_residual {
                return Err(ExactBoundaryRibbonDiagnostic::new(
                    ExactBoundaryRibbonDeclineReason::SupportPcurveStale,
                    context,
                    Some(projected_orientation.max_residual),
                ));
            }
            (&projected_pcurve, projected_orientation)
        }
    };

    let scale = local_sector_width(boundary, other_boundary).min(budget.max_derivative_scale);
    let samples = compute_support_ribbon_samples(SupportRibbonSampleRequest {
        boundary,
        other_boundary,
        surface,
        pcurve,
        parameter: 0.0,
        pcurve_reversed: orientation.reversed,
        derivative_scale: scale,
        budget,
    })
    .map_err(|reason| ExactBoundaryRibbonDiagnostic::new(reason, context, None))?;

    Ok(PreparedExactBoundaryRibbon {
        sector_index: context.sector.sector_index,
        sector_name: context.sector.sector_name.clone(),
        side: context.side,
        edge: context.edge,
        support_half_edge: context.support_half_edge,
        support_face: context.support_face,
        support_family,
        support_surface: surface.clone(),
        support_parameter_reversed: orientation.reversed,
        support_half_edge_curve_forward: context.support_half_edge_curve_forward,
        max_support_residual: orientation.max_residual,
        samples,
    })
}

fn support_half_edge_curve_forward(brep: &BRep, half_edge: HalfEdgeHandle) -> Option<bool> {
    let half_edge_data = brep.half_edges.get(half_edge)?;
    let curve = brep.edges.get(half_edge_data.edge)?.curve.as_ref()?;
    let origin = brep.vertices.get(half_edge_data.origin)?.position;
    let start = curve.geometry.point_at(curve.t_start);
    let end = curve.geometry.point_at(curve.t_end);
    Some(origin.distance_to(&start) <= origin.distance_to(&end))
}

fn boundary_for(
    boundaries: &PreparedExactBoundaryCurves,
    sector_index: usize,
    side: ExactBoundaryOuterSide,
) -> Option<&PreparedExactBoundaryCurve> {
    boundaries
        .curves
        .iter()
        .find(|curve| curve.sector_index == sector_index && curve.side == side)
}

fn other_boundary_for(
    boundaries: &PreparedExactBoundaryCurves,
    sector_index: usize,
    side: ExactBoundaryOuterSide,
) -> Option<&PreparedExactBoundaryCurve> {
    boundaries
        .curves
        .iter()
        .find(|curve| curve.sector_index == sector_index && curve.side != side)
}

fn support_surface(brep: &BRep, context: &RibbonContext<'_>) -> Option<SurfaceGeometry> {
    if let Some(surface) = brep
        .faces
        .get(context.support_face)
        .and_then(|face_data| face_data.surface.as_ref())
        .map(|surface| surface.geometry.clone())
    {
        return Some(surface);
    }
    unique_analytic_surface_for_edge(brep, context)
}

fn unique_analytic_surface_for_edge(
    brep: &BRep,
    context: &RibbonContext<'_>,
) -> Option<SurfaceGeometry> {
    let edge = brep.edges.get(context.edge)?;
    let curve = edge.curve.as_ref()?;
    let tolerance = brep.tolerance.linear().max(1.0e-6);
    let candidates = brep
        .faces
        .iter()
        .filter_map(|(face, face_data)| {
            if face == context.sector.cap_face || face == context.support_face {
                return None;
            }
            let surface = face_data.surface.as_ref()?;
            if matches!(
                surface.geometry,
                SurfaceGeometry::Gregory(_) | SurfaceGeometry::Nurbs(_)
            ) {
                return None;
            }
            let residual = edge_surface_deviation(
                &curve.geometry,
                (curve.t_start, curve.t_end),
                &surface.geometry,
                17,
            )
            .max_deviation;
            (residual <= tolerance).then_some(surface.geometry.clone())
        })
        .collect::<Vec<_>>();
    let mut preferred = candidates
        .iter()
        .filter(|surface| match context.side {
            ExactBoundaryOuterSide::FilletEnd => {
                matches!(surface, SurfaceGeometry::Cylinder { .. })
            }
            ExactBoundaryOuterSide::FaceEdge => matches!(surface, SurfaceGeometry::Plane { .. }),
        })
        .cloned()
        .collect::<Vec<_>>();
    if preferred.len() == 1 {
        return preferred.pop();
    }
    if let Some(surface) = geometrically_unique_surface(&preferred, tolerance) {
        return Some(surface);
    }
    let mut candidates = candidates;
    if candidates.len() == 1 {
        candidates.pop()
    } else {
        geometrically_unique_surface(&candidates, tolerance)
    }
}
