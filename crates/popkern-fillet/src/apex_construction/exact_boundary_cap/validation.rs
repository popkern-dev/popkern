//! Validation gates for exact-boundary apex cap candidates.

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::healing::deviation::edge_surface_deviation;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

use super::corner::PreparedExactBoundaryInnerSpokes;
use super::pcurves::PreparedExactBoundaryPcurves;
use super::ribbon::{PreparedExactBoundaryRibbon, PreparedExactBoundaryRibbons};
use super::sector::ExactBoundaryCapSectorInput;
use super::solve::{ExactBoundarySectorSurface, ExactBoundarySolvedSurfaces};
use super::topology::ExactBoundarySectorTopology;

#[path = "validation_support_normal.rs"]
mod validation_support_normal;
use validation_support_normal::validate_support_normal_continuity;

const DEFAULT_MAX_BOUNDARY_RESIDUAL: f64 = 1.0e-6;
const DEFAULT_MAX_INTERNAL_G0_GAP: f64 = 1.0e-6;
const DEFAULT_MAX_SUPPORT_NORMAL_DEFICIT: f64 = 5.0e-2;
const DEFAULT_MIN_NORMAL_LENGTH: f64 = 1.0e-12;

/// Local validation budgets for exact-boundary candidate acceptance.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct ExactBoundaryValidationBudget {
    /// Maximum cap-side pcurve image residual against locked outer boundaries.
    pub max_outer_boundary_residual: f64,
    /// Maximum support-side pcurve image residual against locked outer boundaries.
    pub max_support_boundary_residual: f64,
    /// Maximum cap/support image disagreement along shared outer seams.
    pub max_pair_boundary_residual: f64,
    /// Maximum gap where generated inner spokes meet the shared pole.
    pub max_internal_g0_gap: f64,
    /// Maximum cap/support unit-normal coplanarity deficit along outer seams.
    pub max_support_normal_deficit: f64,
    /// Minimum sampled normal length for gross foldover rejection.
    pub min_normal_length: f64,
}

impl Default for ExactBoundaryValidationBudget {
    fn default() -> Self {
        Self {
            max_outer_boundary_residual: DEFAULT_MAX_BOUNDARY_RESIDUAL,
            max_support_boundary_residual: DEFAULT_MAX_BOUNDARY_RESIDUAL,
            max_pair_boundary_residual: DEFAULT_MAX_BOUNDARY_RESIDUAL,
            max_internal_g0_gap: DEFAULT_MAX_INTERNAL_G0_GAP,
            max_support_normal_deficit: DEFAULT_MAX_SUPPORT_NORMAL_DEFICIT,
            min_normal_length: DEFAULT_MIN_NORMAL_LENGTH,
        }
    }
}

/// Candidate state consumed by exact-boundary validation.
#[derive(Debug, Clone, Copy)]
pub(crate) struct ExactBoundaryValidationInputs<'a> {
    pub(crate) brep: &'a BRep,
    pub(crate) topology: &'a ExactBoundarySectorTopology,
    pub(crate) sectors: &'a [ExactBoundaryCapSectorInput],
    pub(crate) spokes: &'a PreparedExactBoundaryInnerSpokes,
    pub(crate) ribbons: &'a PreparedExactBoundaryRibbons,
    pub(crate) solved: &'a ExactBoundarySolvedSurfaces,
    pub(crate) pcurves: &'a PreparedExactBoundaryPcurves,
    pub(crate) budget: ExactBoundaryValidationBudget,
}

/// Candidate validation summary for route reports and future UX.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct ExactBoundaryValidationReport {
    /// Number of sector faces that passed topology and surface gates.
    pub validated_sector_count: usize,
    /// Number of paired cap/support pcurves that passed residual gates.
    pub validated_pcurve_pair_count: usize,
    /// Number of generated inner spoke edges that remain topologically live.
    pub validated_spoke_edge_count: usize,
}

/// Stable decline reason for exact-boundary candidate validation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ExactBoundaryValidationDeclineReason {
    /// Materialized topology failed BRep validation.
    TopologyInvalid,
    /// Sector, topology, generated-spoke, or solved-surface inventory disagreed.
    SurfaceInventoryMismatch,
    /// Prepared cap/support pcurve pair inventory disagreed with sectors.
    PcurveInventoryMismatch,
    /// Cap pcurve image did not follow the locked outer boundary.
    OuterBoundaryResidualExceeded,
    /// Support pcurve image did not follow the locked outer boundary.
    SupportBoundaryResidualExceeded,
    /// Cap/support pcurve images did not agree along a shared seam.
    PairBoundaryResidualExceeded,
    /// Generated inner spokes did not meet at the shared pole on the solved surface.
    InternalG0Exceeded,
    /// Sampled sector surface normal collapsed or flipped grossly.
    FoldoverDetected,
    /// Live 3-D edge curve does not lie on the solved cap or support surface.
    LiveEdgeSurfaceResidualExceeded,
    /// The constrained support-normal strip band could not satisfy its postcondition.
    ConstrainedStripBandUnsatisfied,
}

impl ExactBoundaryValidationDeclineReason {
    /// Stable diagnostic code for tests and future UI/API consumers.
    pub(crate) fn code(self) -> &'static str {
        match self {
            Self::TopologyInvalid => "exact_boundary_cap.validation.topology_invalid",
            Self::SurfaceInventoryMismatch => {
                "exact_boundary_cap.validation.surface_inventory_mismatch"
            }
            Self::PcurveInventoryMismatch => {
                "exact_boundary_cap.validation.pcurve_inventory_mismatch"
            }
            Self::OuterBoundaryResidualExceeded => {
                "exact_boundary_cap.validation.outer_boundary_residual_exceeded"
            }
            Self::SupportBoundaryResidualExceeded => {
                "exact_boundary_cap.validation.support_boundary_residual_exceeded"
            }
            Self::PairBoundaryResidualExceeded => {
                "exact_boundary_cap.validation.pair_boundary_residual_exceeded"
            }
            Self::InternalG0Exceeded => "exact_boundary_cap.validation.internal_g0_exceeded",
            Self::FoldoverDetected => "exact_boundary_cap.validation.foldover_detected",
            Self::LiveEdgeSurfaceResidualExceeded => {
                "exact_boundary_cap.validation.live_edge_surface_residual_exceeded"
            }
            Self::ConstrainedStripBandUnsatisfied => {
                "exact_boundary_cap.validation.constrained_strip_band_unsatisfied"
            }
        }
    }
}

/// Structured exact-boundary validation decline.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ExactBoundaryValidationDiagnostic {
    /// Stable decline reason.
    pub reason: ExactBoundaryValidationDeclineReason,
    /// Stable diagnostic code.
    pub code: &'static str,
    pub sector_name: Option<String>,
    /// Measured residual when available.
    pub residual: Option<f64>,
    /// Source edge being validated.
    pub edge: Option<EdgeHandle>,
    /// Support half-edge being validated.
    pub support_half_edge: Option<HalfEdgeHandle>,
    pub support_face: Option<FaceHandle>,
    pub outer_side: Option<String>,
    pub sample_parameter: Option<f64>,
}

impl ExactBoundaryValidationDiagnostic {
    pub(super) fn new(
        reason: ExactBoundaryValidationDeclineReason,
        sector_name: Option<String>,
        residual: Option<f64>,
    ) -> Self {
        Self {
            reason,
            code: reason.code(),
            sector_name,
            residual,
            edge: None,
            support_half_edge: None,
            support_face: None,
            outer_side: None,
            sample_parameter: None,
        }
    }

    fn live_edge_surface_residual(
        sector: &ExactBoundaryCapSectorInput,
        seam: &super::sector::ExactBoundaryOuterSeam,
        residual: f64,
    ) -> Self {
        Self {
            reason: ExactBoundaryValidationDeclineReason::LiveEdgeSurfaceResidualExceeded,
            code: ExactBoundaryValidationDeclineReason::LiveEdgeSurfaceResidualExceeded.code(),
            sector_name: Some(sector.sector_name.clone()),
            residual: Some(residual),
            edge: Some(seam.edge),
            support_half_edge: Some(seam.support_half_edge),
            support_face: Some(seam.support_face),
            outer_side: Some(format!("{:?}", seam.side)),
            sample_parameter: None,
        }
    }

    pub(super) fn constrained_strip_band(
        pair: &super::pcurves::PreparedExactBoundaryPcurvePair,
        residual: f64,
        sample_parameter: f64,
    ) -> Self {
        Self {
            reason: ExactBoundaryValidationDeclineReason::ConstrainedStripBandUnsatisfied,
            code: ExactBoundaryValidationDeclineReason::ConstrainedStripBandUnsatisfied.code(),
            sector_name: Some(pair.sector_name.clone()),
            residual: Some(residual),
            edge: Some(pair.edge),
            support_half_edge: Some(pair.support_half_edge),
            support_face: Some(pair.support_face),
            outer_side: Some(format!("{:?}", pair.side)),
            sample_parameter: Some(sample_parameter),
        }
    }

    pub(crate) fn into_apex_diagnostic(self) -> crate::types::ApexConstructionDiagnostic {
        crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(
            crate::types::ExactBoundaryCapDiagnostic {
                stage: "exact_boundary_cap.validation".to_string(),
                code: self.code.to_string(),
                sector_name: self
                    .sector_name
                    .unwrap_or_else(|| "exact_boundary_cap".to_string()),
                edge: self.edge,
                support_half_edge: self.support_half_edge,
                support_face: self.support_face,
                residual: self.residual,
                outer_side: self.outer_side,
                sample_parameter: self.sample_parameter,
                shared_fact_ids: Vec::new(),
            },
        )
    }
}

/// Validate exact-boundary candidate topology, surfaces, and pcurves before commit.
pub(crate) fn validate_exact_boundary_candidate(
    inputs: ExactBoundaryValidationInputs<'_>,
) -> Result<ExactBoundaryValidationReport, ExactBoundaryValidationDiagnostic> {
    let ExactBoundaryValidationInputs {
        brep,
        topology,
        sectors,
        spokes,
        ribbons,
        solved,
        pcurves,
        budget,
    } = inputs;
    if let Err(err) = brep.validate() {
        return Err(ExactBoundaryValidationDiagnostic::new(
            ExactBoundaryValidationDeclineReason::TopologyInvalid,
            None,
            Some(err.to_string().len() as f64),
        ));
    }
    validate_surface_inventory(topology, sectors, spokes, solved)?;
    validate_pcurve_inventory(sectors, pcurves)?;
    validate_pcurve_residuals(pcurves, budget)?;
    validate_live_edge_surface_residuals(brep, sectors, solved, ribbons, budget)?;
    validate_internal_spokes(spokes, solved, budget)?;
    validate_surface_normals(solved, budget)?;
    validate_support_normal_continuity(solved, pcurves, ribbons, budget)?;
    Ok(ExactBoundaryValidationReport {
        validated_sector_count: sectors.len(),
        validated_pcurve_pair_count: pcurves.pairs.len(),
        validated_spoke_edge_count: topology.spoke_edges.len(),
    })
}

fn validate_live_edge_surface_residuals(
    brep: &BRep,
    sectors: &[ExactBoundaryCapSectorInput],
    solved: &ExactBoundarySolvedSurfaces,
    ribbons: &PreparedExactBoundaryRibbons,
    budget: ExactBoundaryValidationBudget,
) -> Result<(), ExactBoundaryValidationDiagnostic> {
    for sector in sectors {
        let solved_surface = solved_surface(solved, sector).ok_or_else(|| {
            ExactBoundaryValidationDiagnostic::new(
                ExactBoundaryValidationDeclineReason::SurfaceInventoryMismatch,
                Some(sector.sector_name.clone()),
                None,
            )
        })?;
        let cap_surface = SurfaceGeometry::Nurbs(solved_surface.surface.clone());
        for seam in &sector.outer_seams {
            let Some(curve) = brep
                .edges
                .get(seam.edge)
                .and_then(|edge| edge.curve.as_ref())
            else {
                return Err(
                    ExactBoundaryValidationDiagnostic::live_edge_surface_residual(
                        sector,
                        seam,
                        f64::INFINITY,
                    ),
                );
            };
            let borrowed_support_surface;
            let support_surface = if let Some(surface) = brep
                .faces
                .get(seam.support_face)
                .and_then(|face| face.surface.as_ref())
                .map(|surface| &surface.geometry)
            {
                surface
            } else if let Some(ribbon) = ribbon_for_seam(ribbons, sector, seam) {
                let _ribbon_evidence = (
                    &ribbon.sector_name,
                    ribbon.support_parameter_reversed,
                    ribbon.max_support_residual,
                );
                borrowed_support_surface = ribbon.support_surface.clone();
                &borrowed_support_surface
            } else {
                return Err(
                    ExactBoundaryValidationDiagnostic::live_edge_surface_residual(
                        sector,
                        seam,
                        f64::INFINITY,
                    ),
                );
            };
            let t_range = (curve.t_start, curve.t_end);
            let cap_residual =
                edge_surface_deviation(&curve.geometry, t_range, &cap_surface, 9).max_deviation;
            let support_residual =
                edge_surface_deviation(&curve.geometry, t_range, support_surface, 9).max_deviation;
            let residual = cap_residual.max(support_residual);
            if residual > budget.max_pair_boundary_residual {
                return Err(
                    ExactBoundaryValidationDiagnostic::live_edge_surface_residual(
                        sector, seam, residual,
                    ),
                );
            }
        }
    }
    Ok(())
}

fn ribbon_for_seam<'a>(
    ribbons: &'a PreparedExactBoundaryRibbons,
    sector: &ExactBoundaryCapSectorInput,
    seam: &super::sector::ExactBoundaryOuterSeam,
) -> Option<&'a PreparedExactBoundaryRibbon> {
    ribbons.ribbons.iter().find(|ribbon| {
        ribbon.sector_index == sector.sector_index
            && ribbon.side == seam.side
            && ribbon.edge == seam.edge
            && ribbon.support_half_edge == seam.support_half_edge
    })
}

fn validate_surface_inventory(
    topology: &ExactBoundarySectorTopology,
    sectors: &[ExactBoundaryCapSectorInput],
    spokes: &PreparedExactBoundaryInnerSpokes,
    solved: &ExactBoundarySolvedSurfaces,
) -> Result<(), ExactBoundaryValidationDiagnostic> {
    if sectors.is_empty()
        || topology.sectors.len() != sectors.len()
        || topology.spoke_edges.len() != sectors.len()
        || spokes.sectors.len() != sectors.len()
        || solved.validated_surface_count() != sectors.len()
    {
        return Err(ExactBoundaryValidationDiagnostic::new(
            ExactBoundaryValidationDeclineReason::SurfaceInventoryMismatch,
            None,
            None,
        ));
    }
    for sector in sectors {
        if solved_surface(solved, sector).is_none()
            || !topology
                .sectors
                .iter()
                .any(|candidate| candidate.sector_index == sector.sector_index)
        {
            return Err(ExactBoundaryValidationDiagnostic::new(
                ExactBoundaryValidationDeclineReason::SurfaceInventoryMismatch,
                Some(sector.sector_name.clone()),
                None,
            ));
        }
    }
    Ok(())
}

fn validate_pcurve_inventory(
    sectors: &[ExactBoundaryCapSectorInput],
    pcurves: &PreparedExactBoundaryPcurves,
) -> Result<(), ExactBoundaryValidationDiagnostic> {
    if pcurves.validated_pair_count() != sectors.len() * 2 {
        return Err(ExactBoundaryValidationDiagnostic::new(
            ExactBoundaryValidationDeclineReason::PcurveInventoryMismatch,
            None,
            None,
        ));
    }
    for sector in sectors {
        if pcurves
            .pairs
            .iter()
            .filter(|pair| pair.sector_index == sector.sector_index)
            .count()
            != 2
        {
            return Err(ExactBoundaryValidationDiagnostic::new(
                ExactBoundaryValidationDeclineReason::PcurveInventoryMismatch,
                Some(sector.sector_name.clone()),
                None,
            ));
        }
    }
    Ok(())
}

fn validate_pcurve_residuals(
    pcurves: &PreparedExactBoundaryPcurves,
    budget: ExactBoundaryValidationBudget,
) -> Result<(), ExactBoundaryValidationDiagnostic> {
    for pair in &pcurves.pairs {
        if pair.max_cap_residual > budget.max_outer_boundary_residual {
            return Err(ExactBoundaryValidationDiagnostic::new(
                ExactBoundaryValidationDeclineReason::OuterBoundaryResidualExceeded,
                Some(pair.sector_name.clone()),
                Some(pair.max_cap_residual),
            ));
        }
        if pair.max_support_residual > budget.max_support_boundary_residual {
            return Err(ExactBoundaryValidationDiagnostic::new(
                ExactBoundaryValidationDeclineReason::SupportBoundaryResidualExceeded,
                Some(pair.sector_name.clone()),
                Some(pair.max_support_residual),
            ));
        }
        if pair.max_pair_residual > budget.max_pair_boundary_residual {
            return Err(ExactBoundaryValidationDiagnostic::new(
                ExactBoundaryValidationDeclineReason::PairBoundaryResidualExceeded,
                Some(pair.sector_name.clone()),
                Some(pair.max_pair_residual),
            ));
        }
    }
    Ok(())
}

fn validate_internal_spokes(
    spokes: &PreparedExactBoundaryInnerSpokes,
    solved: &ExactBoundarySolvedSurfaces,
    budget: ExactBoundaryValidationBudget,
) -> Result<(), ExactBoundaryValidationDiagnostic> {
    for sector in &spokes.sectors {
        let surface = solved
            .surfaces
            .iter()
            .find(|surface| surface.sector_index == sector.sector_index)
            .ok_or_else(|| {
                ExactBoundaryValidationDiagnostic::new(
                    ExactBoundaryValidationDeclineReason::SurfaceInventoryMismatch,
                    Some(sector.sector_name.clone()),
                    None,
                )
            })?;
        let pole_gap = surface.surface.point_at(1.0, 1.0).distance_to(&spokes.pole);
        if pole_gap > budget.max_internal_g0_gap {
            return Err(ExactBoundaryValidationDiagnostic::new(
                ExactBoundaryValidationDeclineReason::InternalG0Exceeded,
                Some(sector.sector_name.clone()),
                Some(pole_gap),
            ));
        }
    }
    Ok(())
}

fn validate_surface_normals(
    solved: &ExactBoundarySolvedSurfaces,
    budget: ExactBoundaryValidationBudget,
) -> Result<(), ExactBoundaryValidationDiagnostic> {
    for surface in &solved.surfaces {
        for (u, v) in [(0.25, 0.25), (0.5, 0.5), (0.75, 0.75)] {
            let normal_length = surface.surface.normal_at(u, v).length();
            if normal_length < budget.min_normal_length || !normal_length.is_finite() {
                return Err(ExactBoundaryValidationDiagnostic::new(
                    ExactBoundaryValidationDeclineReason::FoldoverDetected,
                    Some(surface.sector_name.clone()),
                    Some(normal_length),
                ));
            }
        }
    }
    Ok(())
}

fn solved_surface<'a>(
    solved: &'a ExactBoundarySolvedSurfaces,
    sector: &ExactBoundaryCapSectorInput,
) -> Option<&'a ExactBoundarySectorSurface> {
    solved
        .surfaces
        .iter()
        .find(|surface| surface.sector_index == sector.sector_index)
}
