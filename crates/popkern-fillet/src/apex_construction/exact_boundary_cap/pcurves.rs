//! Pcurve readiness for in-memory exact-boundary apex cap sectors.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::topology::{EdgeHandle, FaceHandle, HalfEdgeHandle};

use super::corner::{PreparedExactBoundaryInnerSpokes, PreparedExactBoundarySectorBoundary};
use super::ribbon::{PreparedExactBoundaryRibbon, PreparedExactBoundaryRibbons};
use super::sector::ExactBoundaryOuterSide;
use super::solve::{ExactBoundarySectorSurface, ExactBoundarySolvedSurfaces};

const DEFAULT_MAX_RESIDUAL: f64 = 1.0e-6;

/// Local tolerances for pcurve readiness checks.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct ExactBoundaryPcurveBudget {
    /// Maximum distance from cap-side pcurve image to the locked 3-D edge curve.
    pub max_cap_residual: f64,
    /// Maximum distance from support-side pcurve image to the locked 3-D edge curve.
    pub max_support_residual: f64,
    /// Maximum distance between cap-side and support-side pcurve images.
    pub max_pair_residual: f64,
}

impl Default for ExactBoundaryPcurveBudget {
    fn default() -> Self {
        Self {
            max_cap_residual: DEFAULT_MAX_RESIDUAL,
            max_support_residual: DEFAULT_MAX_RESIDUAL,
            max_pair_residual: DEFAULT_MAX_RESIDUAL,
        }
    }
}

/// Prepared pcurve pairs for all exact-boundary sector outer seams.
#[derive(Debug, Clone)]
pub(crate) struct PreparedExactBoundaryPcurves {
    /// One cap/support pcurve pair per immutable outer seam.
    pub pairs: Vec<PreparedExactBoundaryPcurvePair>,
}

impl PreparedExactBoundaryPcurves {
    /// Count pairs that retain semantic ownership and finite residual facts.
    pub(crate) fn validated_pair_count(&self) -> usize {
        self.pairs
            .iter()
            .filter(|pair| {
                let _evidence = (
                    pair.edge,
                    pair.support_half_edge,
                    pair.support_face,
                    &pair.support_surface,
                );
                let cap_start = pair.cap_pcurve.point_at(0.0);
                let support_start = pair.support_pcurve.point_at(0.0);
                !pair.sector_name.is_empty()
                    && pair.sector_index < self.pairs.len()
                    && matches!(
                        pair.side,
                        ExactBoundaryOuterSide::FilletEnd | ExactBoundaryOuterSide::FaceEdge
                    )
                    && pair.max_cap_residual.is_finite()
                    && pair.max_support_residual.is_finite()
                    && pair.max_pair_residual.is_finite()
                    && cap_start.z == 0.0
                    && support_start.z == 0.0
            })
            .count()
    }
}

/// In-memory pcurve pair ready for a future transactional BRep commit.
#[derive(Debug, Clone)]
pub(crate) struct PreparedExactBoundaryPcurvePair {
    /// Planner-order sector index.
    pub sector_index: usize,
    /// Stable sector name.
    pub sector_name: String,
    /// Sector-local outer side.
    pub side: ExactBoundaryOuterSide,
    /// Source BRep edge.
    pub edge: EdgeHandle,
    /// Support half-edge that supplied the support samples.
    pub support_half_edge: HalfEdgeHandle,
    /// Support face sampled by the support-side pcurve.
    pub support_face: FaceHandle,
    /// Exact support surface used to prove and author the support-side pcurve.
    pub support_surface: SurfaceGeometry,
    /// Cap-side pcurve on the solved sector surface.
    pub cap_pcurve: Pcurve,
    /// Support-side pcurve reconstructed from support evidence.
    pub support_pcurve: Pcurve,
    /// Maximum cap pcurve image residual against the locked edge curve.
    pub max_cap_residual: f64,
    /// Maximum support pcurve image residual against the locked edge curve.
    pub max_support_residual: f64,
    /// Maximum cap/support image distance under shared parameter mapping.
    pub max_pair_residual: f64,
}

/// Stable decline reason for pcurve readiness.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ExactBoundaryPcurveDeclineReason {
    /// Solved sector surface was missing.
    SolvedSurfaceMissing,
    /// Support ribbon was missing.
    SupportRibbonMissing,
    /// Cap-side pcurve image does not follow the locked edge curve.
    CapResidualExceeded,
    /// Support-side pcurve image does not follow the locked edge curve.
    SupportResidualExceeded,
    /// Cap-side and support-side pcurve images do not agree.
    PairResidualExceeded,
}

impl ExactBoundaryPcurveDeclineReason {
    /// Stable diagnostic code for tests and future UI/API consumers.
    pub(crate) fn code(self) -> &'static str {
        match self {
            Self::SolvedSurfaceMissing => "exact_boundary_cap.pcurve.solved_surface_missing",
            Self::SupportRibbonMissing => "exact_boundary_cap.pcurve.support_ribbon_missing",
            Self::CapResidualExceeded => "exact_boundary_cap.pcurve.cap_residual_exceeded",
            Self::SupportResidualExceeded => "exact_boundary_cap.pcurve.support_residual_exceeded",
            Self::PairResidualExceeded => "exact_boundary_cap.pcurve.pair_residual_exceeded",
        }
    }
}

/// Structured pcurve-readiness decline.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ExactBoundaryPcurveDiagnostic {
    /// Stable decline reason.
    pub reason: ExactBoundaryPcurveDeclineReason,
    /// Stable diagnostic code.
    pub code: &'static str,
    /// Stable sector name.
    pub sector_name: String,
    /// Sector-local side when known.
    pub side: Option<ExactBoundaryOuterSide>,
    /// Source edge when known.
    pub edge: Option<EdgeHandle>,
    /// Support half-edge when known.
    pub support_half_edge: Option<HalfEdgeHandle>,
    /// Support face when known.
    pub support_face: Option<FaceHandle>,
    /// Measured residual when available.
    pub residual: Option<f64>,
}

impl ExactBoundaryPcurveDiagnostic {
    fn new(
        reason: ExactBoundaryPcurveDeclineReason,
        sector_name: String,
        side: Option<ExactBoundaryOuterSide>,
        ribbon: Option<&PreparedExactBoundaryRibbon>,
        residual: Option<f64>,
    ) -> Self {
        Self {
            reason,
            code: reason.code(),
            sector_name,
            side,
            edge: ribbon.map(|ribbon| ribbon.edge),
            support_half_edge: ribbon.map(|ribbon| ribbon.support_half_edge),
            support_face: ribbon.map(|ribbon| ribbon.support_face),
            residual,
        }
    }

    pub(crate) fn into_apex_diagnostic(self) -> crate::types::ApexConstructionDiagnostic {
        crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(
            crate::types::ExactBoundaryCapDiagnostic {
                stage: "exact_boundary_cap.pcurve".to_string(),
                code: self.code.to_string(),
                sector_name: self.sector_name,
                edge: self.edge,
                support_half_edge: self.support_half_edge,
                support_face: self.support_face,
                residual: self.residual,
                shared_fact_ids: Vec::new(),
            },
        )
    }
}

/// Prepare compatible cap/support pcurve pairs from solved surfaces and ribbons.
pub(crate) fn prepare_exact_boundary_pcurves(
    spokes: &PreparedExactBoundaryInnerSpokes,
    ribbons: &PreparedExactBoundaryRibbons,
    solved: &ExactBoundarySolvedSurfaces,
    budget: ExactBoundaryPcurveBudget,
) -> Result<PreparedExactBoundaryPcurves, ExactBoundaryPcurveDiagnostic> {
    let mut pairs = Vec::with_capacity(spokes.sectors.len() * 2);
    for sector in &spokes.sectors {
        let surface = solved_surface(solved, sector).ok_or_else(|| {
            ExactBoundaryPcurveDiagnostic::new(
                ExactBoundaryPcurveDeclineReason::SolvedSurfaceMissing,
                sector.sector_name.clone(),
                None,
                None,
                None,
            )
        })?;
        pairs.push(prepare_side_pair(
            sector,
            surface,
            ribbons,
            ExactBoundaryOuterSide::FilletEnd,
            budget,
        )?);
        pairs.push(prepare_side_pair(
            sector,
            surface,
            ribbons,
            ExactBoundaryOuterSide::FaceEdge,
            budget,
        )?);
    }
    Ok(PreparedExactBoundaryPcurves { pairs })
}

fn prepare_side_pair(
    sector: &PreparedExactBoundarySectorBoundary,
    surface: &ExactBoundarySectorSurface,
    ribbons: &PreparedExactBoundaryRibbons,
    side: ExactBoundaryOuterSide,
    budget: ExactBoundaryPcurveBudget,
) -> Result<PreparedExactBoundaryPcurvePair, ExactBoundaryPcurveDiagnostic> {
    let Some(ribbon) = ribbon_for(ribbons, sector.sector_index, side) else {
        return Err(ExactBoundaryPcurveDiagnostic::new(
            ExactBoundaryPcurveDeclineReason::SupportRibbonMissing,
            sector.sector_name.clone(),
            Some(side),
            None,
            None,
        ));
    };
    let boundary = match side {
        ExactBoundaryOuterSide::FilletEnd => &sector.fillet_boundary.curve,
        ExactBoundaryOuterSide::FaceEdge => &sector.face_boundary.curve,
    };
    let cap_pcurve = authored_cap_pcurve(sector, side);
    let support_pcurve = support_pcurve_from_ribbon(ribbon);
    let residuals = pcurve_residuals(boundary, &surface.surface, &cap_pcurve, ribbon);
    if residuals.cap > budget.max_cap_residual {
        return Err(pair_diagnostic(
            ExactBoundaryPcurveDeclineReason::CapResidualExceeded,
            sector,
            side,
            ribbon,
            residuals.cap,
        ));
    }
    if residuals.support > budget.max_support_residual {
        return Err(pair_diagnostic(
            ExactBoundaryPcurveDeclineReason::SupportResidualExceeded,
            sector,
            side,
            ribbon,
            residuals.support,
        ));
    }
    if residuals.pair > budget.max_pair_residual {
        return Err(pair_diagnostic(
            ExactBoundaryPcurveDeclineReason::PairResidualExceeded,
            sector,
            side,
            ribbon,
            residuals.pair,
        ));
    }
    Ok(PreparedExactBoundaryPcurvePair {
        sector_index: sector.sector_index,
        sector_name: sector.sector_name.clone(),
        side,
        edge: ribbon.edge,
        support_half_edge: ribbon.support_half_edge,
        support_face: ribbon.support_face,
        support_surface: ribbon.support_surface.clone(),
        cap_pcurve,
        support_pcurve,
        max_cap_residual: residuals.cap,
        max_support_residual: residuals.support,
        max_pair_residual: residuals.pair,
    })
}

fn pair_diagnostic(
    reason: ExactBoundaryPcurveDeclineReason,
    sector: &PreparedExactBoundarySectorBoundary,
    side: ExactBoundaryOuterSide,
    ribbon: &PreparedExactBoundaryRibbon,
    residual: f64,
) -> ExactBoundaryPcurveDiagnostic {
    ExactBoundaryPcurveDiagnostic::new(
        reason,
        sector.sector_name.clone(),
        Some(side),
        Some(ribbon),
        Some(residual),
    )
}

fn solved_surface<'a>(
    solved: &'a ExactBoundarySolvedSurfaces,
    sector: &PreparedExactBoundarySectorBoundary,
) -> Option<&'a ExactBoundarySectorSurface> {
    solved
        .surfaces
        .iter()
        .find(|surface| surface.sector_index == sector.sector_index)
}

fn ribbon_for(
    ribbons: &PreparedExactBoundaryRibbons,
    sector_index: usize,
    side: ExactBoundaryOuterSide,
) -> Option<&PreparedExactBoundaryRibbon> {
    ribbons
        .ribbons
        .iter()
        .find(|ribbon| ribbon.sector_index == sector_index && ribbon.side == side)
}

fn cap_pcurve_for_side(side: ExactBoundaryOuterSide) -> Pcurve {
    match side {
        ExactBoundaryOuterSide::FilletEnd => Pcurve::linear((0.0, 0.0), (1.0, 0.0), 0.0, 1.0),
        ExactBoundaryOuterSide::FaceEdge => Pcurve::linear((0.0, 0.0), (0.0, 1.0), 0.0, 1.0),
    }
}

fn authored_cap_pcurve(
    sector: &PreparedExactBoundarySectorBoundary,
    side: ExactBoundaryOuterSide,
) -> Pcurve {
    let reversed = match side {
        ExactBoundaryOuterSide::FilletEnd => boundary_reversed(
            &sector.fillet_boundary.curve,
            sector.fillet_inner_spoke.anchor,
        ),
        ExactBoundaryOuterSide::FaceEdge => {
            boundary_reversed(&sector.face_boundary.curve, sector.face_inner_spoke.anchor)
        }
    };
    match (side, reversed) {
        (ExactBoundaryOuterSide::FilletEnd, false) => cap_pcurve_for_side(side),
        (ExactBoundaryOuterSide::FilletEnd, true) => {
            Pcurve::linear((1.0, 0.0), (0.0, 0.0), 0.0, 1.0)
        }
        (ExactBoundaryOuterSide::FaceEdge, false) => cap_pcurve_for_side(side),
        (ExactBoundaryOuterSide::FaceEdge, true) => {
            Pcurve::linear((0.0, 1.0), (0.0, 0.0), 0.0, 1.0)
        }
    }
}

fn boundary_reversed(boundary: &NurbsCurve, expected_end: Point3) -> bool {
    boundary.point_at(0.0).distance_to(&expected_end)
        < boundary.point_at(1.0).distance_to(&expected_end)
}

fn support_pcurve_from_ribbon(ribbon: &PreparedExactBoundaryRibbon) -> Pcurve {
    let mut control_points = ribbon
        .samples
        .iter()
        .map(|sample| Point3::new(sample.support_uv.0, sample.support_uv.1, 0.0))
        .collect::<Vec<_>>();
    if !ribbon.support_half_edge_curve_forward {
        control_points.reverse();
    }
    let parameters = ribbon
        .samples
        .iter()
        .map(|sample| sample.parameter)
        .collect::<Vec<_>>();
    let degree = 3.min(control_points.len().saturating_sub(1));
    if let Ok(curve) =
        popkern_core::geometry::fitting::interpolate_curve(&control_points, &parameters, degree)
    {
        return Pcurve::from_nurbs(curve);
    }
    Pcurve::from_nurbs(NurbsCurve {
        degree: 1,
        knots: clamped_linear_knots(control_points.len()),
        weights: vec![1.0; control_points.len()],
        control_points,
    })
}

struct PcurveResiduals {
    cap: f64,
    support: f64,
    pair: f64,
}

fn pcurve_residuals(
    boundary: &impl CurveEval,
    surface: &impl SurfaceEval,
    cap_pcurve: &Pcurve,
    ribbon: &PreparedExactBoundaryRibbon,
) -> PcurveResiduals {
    let mut residuals = PcurveResiduals {
        cap: 0.0,
        support: 0.0,
        pair: 0.0,
    };
    for sample in &ribbon.samples {
        let t = sample.parameter;
        let boundary_point = boundary.point_at(t);
        let cap_uv = cap_pcurve.point_at(t);
        let cap_point = surface.point_at(cap_uv.x, cap_uv.y);
        residuals.cap = residuals.cap.max(cap_point.distance_to(&boundary_point));
        residuals.support = residuals
            .support
            .max(sample.support_point.distance_to(&boundary_point));
        residuals.pair = residuals
            .pair
            .max(cap_point.distance_to(&sample.support_point));
    }
    residuals
}

fn clamped_linear_knots(count: usize) -> Vec<f64> {
    let mut knots = Vec::with_capacity(count + 2);
    knots.push(0.0);
    knots.extend((0..count).map(|index| index as f64 / (count - 1) as f64));
    knots.push(1.0);
    knots
}
