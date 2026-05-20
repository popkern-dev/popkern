//! Support derivative ribbon preparation for exact-boundary apex caps.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

use super::boundary::{PreparedExactBoundaryCurve, PreparedExactBoundaryCurves};
use super::sector::{ExactBoundaryCapSectorInput, ExactBoundaryOuterSide};

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
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct PreparedExactBoundaryRibbons {
    /// One ribbon per immutable outer seam.
    pub ribbons: Vec<PreparedExactBoundaryRibbon>,
}

/// Support derivative ribbon for one immutable outer seam.
#[derive(Debug, Clone, PartialEq)]
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
    /// True when support pcurve orientation is opposite the prepared cap curve.
    pub support_parameter_reversed: bool,
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
                edge: self.edge,
                support_half_edge: Some(self.support_half_edge),
                support_face: Some(self.support_face),
                residual: self.residual,
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
            let surface = support_surface(brep, seam.support_face).ok_or_else(|| {
                ExactBoundaryRibbonDiagnostic::new(
                    ExactBoundaryRibbonDeclineReason::SupportSurfaceMissing,
                    &context,
                    None,
                )
            })?;
            let pcurve = brep
                .half_edges
                .get(seam.support_half_edge)
                .and_then(|half_edge| half_edge.pcurve.as_ref())
                .ok_or_else(|| {
                    ExactBoundaryRibbonDiagnostic::new(
                        ExactBoundaryRibbonDeclineReason::SupportPcurveMissing,
                        &context,
                        None,
                    )
                })?;

            ribbons.push(prepare_ribbon_for_seam(
                &context,
                boundary,
                other_boundary,
                surface,
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
    pcurve: &Pcurve,
    budget: ExactBoundaryRibbonBudget,
) -> Result<PreparedExactBoundaryRibbon, ExactBoundaryRibbonDiagnostic> {
    let orientation = choose_support_pcurve_orientation(boundary, surface, pcurve, budget);
    if orientation.max_residual > budget.max_support_residual {
        return Err(ExactBoundaryRibbonDiagnostic::new(
            ExactBoundaryRibbonDeclineReason::SupportPcurveStale,
            context,
            Some(orientation.max_residual),
        ));
    }

    let scale = local_sector_width(boundary, other_boundary).min(budget.max_derivative_scale);
    let mut samples = Vec::with_capacity(sample_count(budget));
    for parameter in sample_parameters(budget) {
        samples.push(
            compute_support_ribbon_sample(
                boundary,
                other_boundary,
                surface,
                pcurve,
                parameter,
                orientation.reversed,
                scale,
                budget,
            )
            .map_err(|reason| ExactBoundaryRibbonDiagnostic::new(reason, context, None))?,
        );
    }

    Ok(PreparedExactBoundaryRibbon {
        sector_index: context.sector.sector_index,
        sector_name: context.sector.sector_name.clone(),
        side: context.side,
        edge: context.edge,
        support_half_edge: context.support_half_edge,
        support_face: context.support_face,
        support_parameter_reversed: orientation.reversed,
        max_support_residual: orientation.max_residual,
        samples,
    })
}

struct SupportPcurveOrientation {
    reversed: bool,
    max_residual: f64,
}

fn choose_support_pcurve_orientation(
    boundary: &PreparedExactBoundaryCurve,
    surface: &SurfaceGeometry,
    pcurve: &Pcurve,
    budget: ExactBoundaryRibbonBudget,
) -> SupportPcurveOrientation {
    let forward = max_support_residual(boundary, surface, pcurve, false, budget);
    let reversed = max_support_residual(boundary, surface, pcurve, true, budget);
    if reversed <= forward {
        SupportPcurveOrientation {
            reversed: true,
            max_residual: reversed,
        }
    } else {
        SupportPcurveOrientation {
            reversed: false,
            max_residual: forward,
        }
    }
}

fn max_support_residual(
    boundary: &PreparedExactBoundaryCurve,
    surface: &SurfaceGeometry,
    pcurve: &Pcurve,
    reversed: bool,
    budget: ExactBoundaryRibbonBudget,
) -> f64 {
    sample_parameters(budget)
        .map(|parameter| {
            let boundary_point = boundary.curve.point_at(parameter);
            let uv = pcurve.point_at(support_pcurve_parameter(pcurve, parameter, reversed));
            let support_point = surface.point_at(uv.x, uv.y);
            boundary_point.distance_to(&support_point)
        })
        .fold(0.0, f64::max)
}

pub(super) fn compute_support_ribbon_sample(
    boundary: &PreparedExactBoundaryCurve,
    other_boundary: &PreparedExactBoundaryCurve,
    surface: &SurfaceGeometry,
    pcurve: &Pcurve,
    parameter: f64,
    pcurve_reversed: bool,
    derivative_scale: f64,
    budget: ExactBoundaryRibbonBudget,
) -> Result<ExactBoundaryRibbonSample, ExactBoundaryRibbonDeclineReason> {
    let pcurve_parameter = support_pcurve_parameter(pcurve, parameter, pcurve_reversed);
    let uv = pcurve.point_at(pcurve_parameter);
    let support_point = surface.point_at(uv.x, uv.y);
    let boundary_point = boundary.curve.point_at(parameter);
    let normal = surface.normal_at(uv.x, uv.y);
    let Some(support_normal) = normal.normalized() else {
        return Err(ExactBoundaryRibbonDeclineReason::DegenerateSupportNormal);
    };
    if normal.length() < budget.min_normal_length {
        return Err(ExactBoundaryRibbonDeclineReason::DegenerateSupportNormal);
    }

    let tangent = boundary.curve.derivative_at(parameter, 1);
    let tangent_in_plane = tangent - support_normal * tangent.dot(&support_normal);
    if tangent_in_plane.length() < budget.min_tangent_length {
        return Err(ExactBoundaryRibbonDeclineReason::TangentAmbiguity);
    }
    let Some(boundary_tangent) = tangent_in_plane.normalized() else {
        return Err(ExactBoundaryRibbonDeclineReason::TangentAmbiguity);
    };

    let cross = support_normal.cross(&boundary_tangent);
    if cross.length() < budget.min_cross_length {
        return Err(ExactBoundaryRibbonDeclineReason::CrossDerivativeDegenerate);
    }
    let Some(mut cross_unit) = cross.normalized() else {
        return Err(ExactBoundaryRibbonDeclineReason::CrossDerivativeDegenerate);
    };
    let inward = other_boundary.curve.point_at(0.5) - boundary_point;
    if cross_unit.dot(&inward) < 0.0 {
        cross_unit = cross_unit * -1.0;
    }
    let cross_derivative = cross_unit * derivative_scale;

    Ok(ExactBoundaryRibbonSample {
        parameter,
        support_uv: (uv.x, uv.y),
        boundary_point,
        support_point,
        boundary_tangent,
        support_normal,
        cross_derivative,
    })
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

fn support_surface(brep: &BRep, face: FaceHandle) -> Option<&SurfaceGeometry> {
    brep.faces
        .get(face)
        .and_then(|face_data| face_data.surface.as_ref())
        .map(|surface| &surface.geometry)
}

fn support_pcurve_parameter(pcurve: &Pcurve, parameter: f64, reversed: bool) -> f64 {
    let mapped = if reversed { 1.0 - parameter } else { parameter };
    pcurve.t_start + (pcurve.t_end - pcurve.t_start) * mapped
}

fn local_sector_width(
    boundary: &PreparedExactBoundaryCurve,
    other_boundary: &PreparedExactBoundaryCurve,
) -> f64 {
    boundary
        .curve
        .point_at(0.5)
        .distance_to(&other_boundary.curve.point_at(0.5))
}

fn sample_count(budget: ExactBoundaryRibbonBudget) -> usize {
    budget.sample_count.max(2)
}

fn sample_parameters(
    budget: ExactBoundaryRibbonBudget,
) -> impl Iterator<Item = f64> + Clone + ExactSizeIterator {
    let count = sample_count(budget);
    (0..count).map(move |index| index as f64 / (count - 1) as f64)
}
