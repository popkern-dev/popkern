//! Generated inner-spoke boundary model for exact-boundary apex cap sectors.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::{NurbsCurve, Point3};

use super::boundary::{PreparedExactBoundaryCurve, PreparedExactBoundaryCurves};
use super::sector::{ExactBoundaryCapSectorInput, ExactBoundaryOuterSide};

const DEFAULT_G0_TOLERANCE: f64 = 1.0e-7;
const DEFAULT_MIN_SPOKE_LENGTH: f64 = 1.0e-6;
const GENERATED_ANCHOR_RECONCILIATION_TOLERANCE: f64 = 5.0e-3;

/// Local tolerances for generated inner-spoke preparation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct ExactBoundaryInnerSpokeBudget {
    /// Maximum distance for identifying shared outer corners and adjacent spoke anchors.
    pub g0_tolerance: f64,
    /// Minimum generated spoke length.
    pub min_spoke_length: f64,
}

impl Default for ExactBoundaryInnerSpokeBudget {
    fn default() -> Self {
        Self {
            g0_tolerance: DEFAULT_G0_TOLERANCE,
            min_spoke_length: DEFAULT_MIN_SPOKE_LENGTH,
        }
    }
}

/// Generated inner-spoke model for all sectors.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct PreparedExactBoundaryInnerSpokes {
    /// One generated sector-boundary record per exact-boundary sector.
    pub sectors: Vec<PreparedExactBoundarySectorBoundary>,
    /// Deterministic pole used by all generated spokes.
    pub pole: Point3,
    /// Blend-ward pull direction for the PK-255 pole search (the
    /// negated sum of the confident end-section plane normals), and
    /// the search scale (the largest fillet-end arc chord). `None`
    /// when no plane orientation was confident.
    pub pole_search: Option<(popkern_core::geometry::Vec3, f64)>,
}

/// Four-sided sector boundary data before surface solve.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct PreparedExactBoundarySectorBoundary {
    /// Planner-order sector index.
    pub sector_index: usize,
    /// Stable sector name.
    pub sector_name: String,
    /// Immutable fillet-end outer boundary.
    pub fillet_boundary: PreparedExactBoundaryCurve,
    /// Immutable face-edge outer boundary. For a welded
    /// shared-corner sector (PK-219) this is the SYNTHESIZED
    /// point-contact boundary at the weld — the PK-125 polar-wedge
    /// degenerate-outer-side machinery drives the solve; commit
    /// skips the edge-less side via the welded seam records.
    pub face_boundary: PreparedExactBoundaryCurve,
    /// Generated spoke from the fillet-side inner anchor to the pole.
    pub fillet_inner_spoke: GeneratedInnerSpoke,
    /// Generated spoke from the face-side inner anchor to the pole.
    pub face_inner_spoke: GeneratedInnerSpoke,
    /// Measured G0 gap where the two immutable outer boundaries meet.
    pub outer_corner_gap: f64,
}

/// Generated cap-internal boundary with provenance.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct GeneratedInnerSpoke {
    /// Stable generated spoke name.
    pub name: String,
    /// Sector-local side that supplied the outer anchor.
    pub source_side: ExactBoundaryOuterSide,
    /// Anchor on the immutable outer boundary.
    pub anchor: Point3,
    /// Shared deterministic pole.
    pub pole: Point3,
    /// Degree-1 generated curve from anchor to pole.
    pub curve: NurbsCurve,
}

/// Stable decline reason for inner-spoke preparation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ExactBoundaryInnerSpokeDeclineReason {
    /// An expected prepared boundary was missing.
    BoundaryMissing,
    /// The two immutable outer boundaries did not share a corner.
    OuterCornerMismatch,
    /// A generated spoke collapsed below the local tolerance.
    CollapsedSpoke,
    /// Adjacent sector spokes did not share the same outer anchor.
    AdjacentSpokeMismatch,
}

impl ExactBoundaryInnerSpokeDeclineReason {
    /// Stable diagnostic code for tests and future UI/API consumers.
    pub(crate) fn code(self) -> &'static str {
        match self {
            Self::BoundaryMissing => "exact_boundary_cap.inner_spoke.boundary_missing",
            Self::OuterCornerMismatch => "exact_boundary_cap.inner_spoke.outer_corner_mismatch",
            Self::CollapsedSpoke => "exact_boundary_cap.inner_spoke.collapsed",
            Self::AdjacentSpokeMismatch => "exact_boundary_cap.inner_spoke.adjacent_g0_mismatch",
        }
    }
}

/// Structured inner-spoke decline.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ExactBoundaryInnerSpokeDiagnostic {
    /// Stable decline reason.
    pub reason: ExactBoundaryInnerSpokeDeclineReason,
    /// Stable diagnostic code.
    pub code: &'static str,
    /// Planner-order sector index when available.
    pub sector_index: Option<usize>,
    /// Stable sector name when available.
    pub sector_name: Option<String>,
    /// Measured residual when available.
    pub residual: Option<f64>,
}

impl ExactBoundaryInnerSpokeDiagnostic {
    fn new(
        reason: ExactBoundaryInnerSpokeDeclineReason,
        sector: Option<&ExactBoundaryCapSectorInput>,
        residual: Option<f64>,
    ) -> Self {
        Self {
            reason,
            code: reason.code(),
            sector_index: sector.map(|sector| sector.sector_index),
            sector_name: sector.map(|sector| sector.sector_name.clone()),
            residual,
        }
    }

    pub(crate) fn into_apex_diagnostic(self) -> crate::types::ApexConstructionDiagnostic {
        crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(
            crate::types::ExactBoundaryCapDiagnostic {
                stage: "exact_boundary_cap.inner_spoke".to_string(),
                code: self.code.to_string(),
                detail_code: None,
                sector_name: self.sector_name.unwrap_or_else(|| "unknown".to_string()),
                edge: None,
                support_half_edge: None,
                support_face: None,
                residual: self.residual,
                outer_side: None,
                sample_parameter: None,
                tolerance: None,
                shared_fact_ids: Vec::new(),
            },
        )
    }
}

#[path = "corner_spoke_shaping.rs"]
mod corner_spoke_shaping;
use corner_spoke_shaping::shape_spokes_for_corner_compatibility;

/// Prepare generated inner spokes and complete four-sided sector boundaries.
/// Prepare the inner spokes; `pole_override` supplies an explicit
/// pole (PK-255 pole-search candidates), `None` keeps the
/// shared-corner centroid.
pub(crate) fn prepare_exact_boundary_inner_spokes_with_pole(
    sectors: &[ExactBoundaryCapSectorInput],
    boundaries: &PreparedExactBoundaryCurves,
    budget: ExactBoundaryInnerSpokeBudget,
    pole_override: Option<Point3>,
) -> Result<PreparedExactBoundaryInnerSpokes, ExactBoundaryInnerSpokeDiagnostic> {
    let mut sector_inputs = Vec::with_capacity(sectors.len());
    let mut shared_corners = Vec::with_capacity(sectors.len());
    for sector in sectors {
        let fillet = boundary_for(
            boundaries,
            sector.sector_index,
            ExactBoundaryOuterSide::FilletEnd,
        )
        .ok_or_else(|| {
            ExactBoundaryInnerSpokeDiagnostic::new(
                ExactBoundaryInnerSpokeDeclineReason::BoundaryMissing,
                Some(sector),
                None,
            )
        })?;
        // Welded sectors (PK-219) carry a SYNTHESIZED point-contact
        // face boundary from curve preparation; the normal
        // classification handles it (shared point, far point, and
        // spoke anchor all resolve to the weld).
        let face = boundary_for(
            boundaries,
            sector.sector_index,
            ExactBoundaryOuterSide::FaceEdge,
        )
        .ok_or_else(|| {
            ExactBoundaryInnerSpokeDiagnostic::new(
                ExactBoundaryInnerSpokeDeclineReason::BoundaryMissing,
                Some(sector),
                None,
            )
        })?;
        let corner = classify_outer_corner(fillet, face);
        // Corner acceptance runs at the ANCHOR RECONCILIATION scale,
        // not the post-reconciliation spoke gate: `shared_point` is
        // already the midpoint of the two curve ends, built exactly
        // to absorb small corner gaps, and the boundary-conditioned
        // fillet marches (chi-calibrated m-smoothing) legitimately
        // move band ends by ~5e-4 (measured: pyramid4_uniform sector
        // 2 gap 5.65e-4 declined the whole exact-boundary route at
        // the old 1e-7 gate and the fallback cap failed G1 at 0.86).
        // Downstream cap acceptance still measures the committed
        // seams at its own commercial gates.
        if corner.gap > GENERATED_ANCHOR_RECONCILIATION_TOLERANCE.max(budget.g0_tolerance) {
            return Err(ExactBoundaryInnerSpokeDiagnostic::new(
                ExactBoundaryInnerSpokeDeclineReason::OuterCornerMismatch,
                Some(sector),
                Some(corner.gap),
            ));
        }
        shared_corners.push(corner.shared_point);
        sector_inputs.push((sector, fillet.clone(), face.clone(), corner));
    }

    let pole = match pole_override {
        Some(pole) => pole,
        None => centroid(&shared_corners).ok_or_else(|| {
            ExactBoundaryInnerSpokeDiagnostic::new(
                ExactBoundaryInnerSpokeDeclineReason::BoundaryMissing,
                None,
                None,
            )
        })?,
    };
    let pole_search = super::corner_pole_search_hints::pole_search_hints(
        &sector_inputs
            .iter()
            .map(|(_, fillet, _, _)| &fillet.curve)
            .collect::<Vec<_>>(),
    );

    let mut prepared = Vec::with_capacity(sector_inputs.len());
    for (sector, fillet, face, corner) in sector_inputs {
        let fillet_inner_spoke = generated_spoke(
            sector,
            ExactBoundaryOuterSide::FilletEnd,
            corner.fillet_far_point,
            pole,
            budget,
        )?;
        let face_inner_spoke = generated_spoke(
            sector,
            ExactBoundaryOuterSide::FaceEdge,
            corner.face_far_point,
            pole,
            budget,
        )?;
        prepared.push(PreparedExactBoundarySectorBoundary {
            sector_index: sector.sector_index,
            sector_name: sector.sector_name.clone(),
            fillet_boundary: fillet,
            face_boundary: face,
            fillet_inner_spoke,
            face_inner_spoke,
            outer_corner_gap: corner.gap,
        });
    }

    reconcile_adjacent_spoke_anchors(&mut prepared, GENERATED_ANCHOR_RECONCILIATION_TOLERANCE);
    shape_spokes_for_corner_compatibility(&mut prepared);
    validate_adjacent_spoke_g0(&prepared, budget)?;
    Ok(PreparedExactBoundaryInnerSpokes {
        sectors: prepared,
        pole,
        pole_search,
    })
}

struct OuterCornerClassification {
    shared_point: Point3,
    fillet_far_point: Point3,
    face_far_point: Point3,
    gap: f64,
}

fn classify_outer_corner(
    fillet: &PreparedExactBoundaryCurve,
    face: &PreparedExactBoundaryCurve,
) -> OuterCornerClassification {
    let fillet_start = fillet.curve.point_at(0.0);
    let fillet_end = fillet.curve.point_at(1.0);
    let face_start = face.curve.point_at(0.0);
    let face_end = face.curve.point_at(1.0);
    let candidates = [
        (fillet_start, fillet_end, face_start, face_end),
        (fillet_start, fillet_end, face_end, face_start),
        (fillet_end, fillet_start, face_start, face_end),
        (fillet_end, fillet_start, face_end, face_start),
    ];
    let mut best = candidates[0];
    let mut best_gap = best.0.distance_to(&best.2);
    for candidate in candidates.iter().copied().skip(1) {
        let gap = candidate.0.distance_to(&candidate.2);
        if gap < best_gap {
            best = candidate;
            best_gap = gap;
        }
    }
    OuterCornerClassification {
        shared_point: midpoint(best.0, best.2),
        fillet_far_point: best.1,
        face_far_point: best.3,
        gap: best_gap,
    }
}

fn generated_spoke(
    sector: &ExactBoundaryCapSectorInput,
    source_side: ExactBoundaryOuterSide,
    anchor: Point3,
    pole: Point3,
    budget: ExactBoundaryInnerSpokeBudget,
) -> Result<GeneratedInnerSpoke, ExactBoundaryInnerSpokeDiagnostic> {
    let length = anchor.distance_to(&pole);
    if length < budget.min_spoke_length {
        return Err(ExactBoundaryInnerSpokeDiagnostic::new(
            ExactBoundaryInnerSpokeDeclineReason::CollapsedSpoke,
            Some(sector),
            Some(length),
        ));
    }
    Ok(GeneratedInnerSpoke {
        name: format!(
            "{}.inner_spoke.{}",
            sector.sector_name,
            match source_side {
                ExactBoundaryOuterSide::FilletEnd => "fillet_anchor_to_pole",
                ExactBoundaryOuterSide::FaceEdge => "face_anchor_to_pole",
            }
        ),
        source_side,
        anchor,
        pole,
        curve: NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![anchor, pole],
            weights: vec![1.0, 1.0],
        },
    })
}

fn validate_adjacent_spoke_g0(
    sectors: &[PreparedExactBoundarySectorBoundary],
    budget: ExactBoundaryInnerSpokeBudget,
) -> Result<(), ExactBoundaryInnerSpokeDiagnostic> {
    if sectors.len() < 2 {
        return Ok(());
    }
    for index in 0..sectors.len() {
        let current = &sectors[index];
        let gap = sectors
            .iter()
            .enumerate()
            .filter(|(candidate_index, _)| *candidate_index != index)
            .map(|(_, candidate)| {
                current
                    .face_inner_spoke
                    .anchor
                    .distance_to(&candidate.fillet_inner_spoke.anchor)
            })
            .fold(f64::INFINITY, f64::min);
        if gap > budget.g0_tolerance {
            return Err(ExactBoundaryInnerSpokeDiagnostic {
                reason: ExactBoundaryInnerSpokeDeclineReason::AdjacentSpokeMismatch,
                code: ExactBoundaryInnerSpokeDeclineReason::AdjacentSpokeMismatch.code(),
                sector_index: Some(current.sector_index),
                sector_name: Some(current.sector_name.clone()),
                residual: Some(gap),
            });
        }
    }
    Ok(())
}

fn reconcile_adjacent_spoke_anchors(
    sectors: &mut [PreparedExactBoundarySectorBoundary],
    tolerance: f64,
) {
    if sectors.len() < 2 {
        return;
    }
    for index in 0..sectors.len() {
        let Some(nearest_index) = nearest_fillet_anchor_index(sectors, index) else {
            continue;
        };
        let face_anchor = sectors[index].face_inner_spoke.anchor;
        let fillet_anchor = sectors[nearest_index].fillet_inner_spoke.anchor;
        if face_anchor.distance_to(&fillet_anchor) > tolerance {
            continue;
        }
        reset_generated_spoke_anchor(&mut sectors[nearest_index].fillet_inner_spoke, face_anchor);
    }
}

fn nearest_fillet_anchor_index(
    sectors: &[PreparedExactBoundarySectorBoundary],
    index: usize,
) -> Option<usize> {
    sectors
        .iter()
        .enumerate()
        .filter(|(candidate_index, _)| *candidate_index != index)
        .min_by(|(_, a), (_, b)| {
            let anchor = sectors[index].face_inner_spoke.anchor;
            anchor
                .distance_to(&a.fillet_inner_spoke.anchor)
                .partial_cmp(&anchor.distance_to(&b.fillet_inner_spoke.anchor))
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(candidate_index, _)| candidate_index)
}

fn reset_generated_spoke_anchor(spoke: &mut GeneratedInnerSpoke, anchor: Point3) {
    spoke.anchor = anchor;
    spoke.curve.control_points = vec![anchor, spoke.pole];
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

fn centroid(points: &[Point3]) -> Option<Point3> {
    if points.is_empty() {
        return None;
    }
    let sum = points.iter().fold((0.0, 0.0, 0.0), |(x, y, z), point| {
        (x + point.x, y + point.y, z + point.z)
    });
    let scale = 1.0 / points.len() as f64;
    Some(Point3::new(sum.0 * scale, sum.1 * scale, sum.2 * scale))
}

fn midpoint(a: Point3, b: Point3) -> Point3 {
    Point3::new((a.x + b.x) * 0.5, (a.y + b.y) * 0.5, (a.z + b.z) * 0.5)
}
