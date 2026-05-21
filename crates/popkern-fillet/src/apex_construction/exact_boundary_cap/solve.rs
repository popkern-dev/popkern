//! In-memory Hermite/Coons solve for exact-boundary apex cap sectors.

use popkern_core::geometry::SurfaceEval;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::{NurbsCurve, NurbsSurface, Point3, Vec3};

use super::corner::{PreparedExactBoundaryInnerSpokes, PreparedExactBoundarySectorBoundary};
use super::ribbon::{PreparedExactBoundaryRibbon, PreparedExactBoundaryRibbons};
use super::sector::ExactBoundaryOuterSide;

const DEFAULT_CORNER_TOLERANCE: f64 = 1.0e-7;
const DEFAULT_MIN_NORMAL_LENGTH: f64 = 1.0e-12;

/// Local tolerances for the in-memory surface solve.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct ExactBoundarySolveBudget {
    /// Maximum corner mismatch between adjacent locked boundaries.
    pub corner_tolerance: f64,
    /// Minimum normal length in gross foldover checks.
    pub min_normal_length: f64,
}

impl Default for ExactBoundarySolveBudget {
    fn default() -> Self {
        Self {
            corner_tolerance: DEFAULT_CORNER_TOLERANCE,
            min_normal_length: DEFAULT_MIN_NORMAL_LENGTH,
        }
    }
}

/// Candidate surfaces built in memory before BRep commit.
#[derive(Debug, Clone)]
pub(crate) struct ExactBoundarySolvedSurfaces {
    /// One candidate surface per sector.
    pub surfaces: Vec<ExactBoundarySectorSurface>,
}

impl ExactBoundarySolvedSurfaces {
    /// Count candidates that retain their semantic identity and rectangular control net.
    pub(crate) fn validated_surface_count(&self) -> usize {
        self.surfaces
            .iter()
            .filter(|surface| {
                !surface.sector_name.is_empty()
                    && surface.surface.control_points.len()
                        == surface.surface.num_u * surface.surface.num_v
                    && surface.max_corner_gap >= 0.0
                    && surface.sector_index < self.surfaces.len()
            })
            .count()
    }
}

/// In-memory exact-boundary sector surface candidate.
#[derive(Debug, Clone)]
pub(crate) struct ExactBoundarySectorSurface {
    /// Planner-order sector index.
    pub sector_index: usize,
    /// Stable sector name.
    pub sector_name: String,
    /// Candidate tensor-product surface.
    pub surface: NurbsSurface,
    /// Maximum locked corner mismatch measured before solving.
    pub max_corner_gap: f64,
}

/// Stable decline reason for the in-memory solve.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ExactBoundarySolveDeclineReason {
    /// Boundary curve compatibility was missing or inconsistent.
    BoundaryCompatibilityMissing,
    /// Locked boundary corners disagree.
    CornerDerivativeConflict,
    /// Gross foldover check found a collapsed or flipped sample.
    GrossFoldover,
}

impl ExactBoundarySolveDeclineReason {
    /// Stable diagnostic code for tests and future UI/API consumers.
    pub(crate) fn code(self) -> &'static str {
        match self {
            Self::BoundaryCompatibilityMissing => {
                "exact_boundary_cap.solve.boundary_compatibility_missing"
            }
            Self::CornerDerivativeConflict => "exact_boundary_cap.corner_derivative_conflict",
            Self::GrossFoldover => "exact_boundary_cap.solve.gross_foldover",
        }
    }
}

/// Structured solve decline.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ExactBoundarySolveDiagnostic {
    /// Stable decline reason.
    pub reason: ExactBoundarySolveDeclineReason,
    /// Stable diagnostic code.
    pub code: &'static str,
    /// Planner-order sector index.
    pub sector_index: usize,
    /// Stable sector name.
    pub sector_name: String,
    /// Measured residual when available.
    pub residual: Option<f64>,
}

impl ExactBoundarySolveDiagnostic {
    fn new(
        reason: ExactBoundarySolveDeclineReason,
        sector: &PreparedExactBoundarySectorBoundary,
        residual: Option<f64>,
    ) -> Self {
        Self {
            reason,
            code: reason.code(),
            sector_index: sector.sector_index,
            sector_name: sector.sector_name.clone(),
            residual,
        }
    }

    pub(crate) fn into_apex_diagnostic(self) -> crate::types::ApexConstructionDiagnostic {
        crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(
            crate::types::ExactBoundaryCapDiagnostic {
                stage: "exact_boundary_cap.solve".to_string(),
                code: self.code.to_string(),
                sector_name: self.sector_name,
                edge: None,
                support_half_edge: None,
                support_face: None,
                residual: self.residual,
                shared_fact_ids: Vec::new(),
            },
        )
    }
}

/// Build deterministic in-memory sector surfaces from locked boundaries.
pub(crate) fn solve_exact_boundary_surfaces(
    spokes: &PreparedExactBoundaryInnerSpokes,
    ribbons: &PreparedExactBoundaryRibbons,
    budget: ExactBoundarySolveBudget,
) -> Result<ExactBoundarySolvedSurfaces, ExactBoundarySolveDiagnostic> {
    let mut surfaces = Vec::with_capacity(spokes.sectors.len());
    for sector in &spokes.sectors {
        surfaces.push(solve_sector(sector, ribbons, budget)?);
    }
    Ok(ExactBoundarySolvedSurfaces { surfaces })
}

fn solve_sector(
    sector: &PreparedExactBoundarySectorBoundary,
    ribbons: &PreparedExactBoundaryRibbons,
    budget: ExactBoundarySolveBudget,
) -> Result<ExactBoundarySectorSurface, ExactBoundarySolveDiagnostic> {
    let fillet = orient_curve(
        &sector.fillet_boundary.curve,
        sector.fillet_inner_spoke.anchor,
    );
    let face = orient_curve(&sector.face_boundary.curve, sector.face_inner_spoke.anchor);
    let Some((num_u, num_v)) = compatible_dimensions(&fillet, &face) else {
        return Err(ExactBoundarySolveDiagnostic::new(
            ExactBoundarySolveDeclineReason::BoundaryCompatibilityMissing,
            sector,
            None,
        ));
    };
    let shared = fillet.point_at(0.0);
    let corner_gap = shared.distance_to(&face.point_at(0.0));
    if corner_gap > budget.corner_tolerance {
        return Err(ExactBoundarySolveDiagnostic::new(
            ExactBoundarySolveDeclineReason::CornerDerivativeConflict,
            sector,
            Some(corner_gap),
        ));
    }

    let top = line_curve_like(&fillet, sector.face_inner_spoke.anchor, sector.pole());
    let right = line_curve_like(&face, sector.fillet_inner_spoke.anchor, sector.pole());
    let mut controls = coons_controls(&fillet, &face, &top, &right, num_u, num_v);
    apply_derivative_rows(sector, ribbons, &mut controls, num_u, num_v);
    let surface = NurbsSurface {
        degree_u: fillet.degree,
        degree_v: face.degree,
        knots_u: fillet.knots.clone(),
        knots_v: face.knots.clone(),
        control_points: controls,
        num_u,
        num_v,
        weights: vec![1.0; num_u * num_v],
    };
    if has_gross_foldover(&surface, budget) {
        return Err(ExactBoundarySolveDiagnostic::new(
            ExactBoundarySolveDeclineReason::GrossFoldover,
            sector,
            None,
        ));
    }
    Ok(ExactBoundarySectorSurface {
        sector_index: sector.sector_index,
        sector_name: sector.sector_name.clone(),
        surface,
        max_corner_gap: corner_gap,
    })
}

trait SectorPole {
    fn pole(&self) -> Point3;
}

impl SectorPole for PreparedExactBoundarySectorBoundary {
    fn pole(&self) -> Point3 {
        self.fillet_inner_spoke.pole
    }
}

fn compatible_dimensions(fillet: &NurbsCurve, face: &NurbsCurve) -> Option<(usize, usize)> {
    if fillet.control_points.len() < 2 || face.control_points.len() < 2 {
        return None;
    }
    Some((fillet.control_points.len(), face.control_points.len()))
}

fn orient_curve(curve: &NurbsCurve, expected_end: Point3) -> NurbsCurve {
    let end_gap = curve.point_at(1.0).distance_to(&expected_end);
    let start_gap = curve.point_at(0.0).distance_to(&expected_end);
    if start_gap < end_gap {
        reverse_curve(curve)
    } else {
        curve.clone()
    }
}

fn line_curve_like(template: &NurbsCurve, start: Point3, end: Point3) -> NurbsCurve {
    let (lo, hi) = template.parameter_range();
    let span = hi - lo;
    let degree = template.degree;
    let cps = (0..template.control_points.len())
        .map(|i| {
            let t = if degree == 0 {
                i as f64 / (template.control_points.len() - 1) as f64
            } else {
                let greville =
                    template.knots[i + 1..=i + degree].iter().sum::<f64>() / degree as f64;
                (greville - lo) / span
            };
            lerp(start, end, t.clamp(0.0, 1.0))
        })
        .collect();
    NurbsCurve {
        degree,
        knots: template.knots.clone(),
        control_points: cps,
        weights: vec![1.0; template.control_points.len()],
    }
}

fn coons_controls(
    bottom: &NurbsCurve,
    left: &NurbsCurve,
    top: &NurbsCurve,
    right: &NurbsCurve,
    num_u: usize,
    num_v: usize,
) -> Vec<Point3> {
    let p00 = bottom.control_points[0];
    let p10 = bottom.control_points[num_u - 1];
    let p01 = left.control_points[num_v - 1];
    let p11 = top.control_points[num_u - 1];
    let mut controls = vec![p00; num_u * num_v];
    for i in 0..num_u {
        let u = i as f64 / (num_u - 1) as f64;
        for j in 0..num_v {
            let v = j as f64 / (num_v - 1) as f64;
            let c_u = blend_point(bottom.control_points[i], top.control_points[i], v);
            let c_v = blend_point(left.control_points[j], right.control_points[j], u);
            let bilinear = bilinear_point(p00, p10, p01, p11, u, v);
            controls[index(i, j, num_v)] =
                add_vec(c_u, c_v) - (bilinear - Point3::new(0.0, 0.0, 0.0));
        }
    }
    controls
}

fn apply_derivative_rows(
    sector: &PreparedExactBoundarySectorBoundary,
    ribbons: &PreparedExactBoundaryRibbons,
    controls: &mut [Point3],
    num_u: usize,
    num_v: usize,
) {
    let mut derivative_controls = vec![None; controls.len()];
    if let Some(ribbon) = ribbon_for(
        ribbons,
        sector.sector_index,
        ExactBoundaryOuterSide::FilletEnd,
    ) {
        for i in 1..num_u.saturating_sub(1) {
            let parameter = i as f64 / (num_u - 1) as f64;
            let derivative = ribbon_derivative(ribbon, parameter);
            let base = controls[index(i, 0, num_v)];
            merge_derivative_control(
                &mut derivative_controls,
                index(i, 1, num_v),
                base + derivative * (1.0 / sector.face_boundary.curve.degree as f64),
            );
        }
    }
    if let Some(ribbon) = ribbon_for(
        ribbons,
        sector.sector_index,
        ExactBoundaryOuterSide::FaceEdge,
    ) {
        for j in 1..num_v.saturating_sub(1) {
            let parameter = j as f64 / (num_v - 1) as f64;
            let derivative = ribbon_derivative(ribbon, parameter);
            let base = controls[index(0, j, num_v)];
            merge_derivative_control(
                &mut derivative_controls,
                index(1, j, num_v),
                base + derivative * (1.0 / sector.fillet_boundary.curve.degree as f64),
            );
        }
    }
    for (control, derivative_control) in controls.iter_mut().zip(derivative_controls) {
        if let Some(derivative_control) = derivative_control {
            *control = derivative_control;
        }
    }
}

fn merge_derivative_control(controls: &mut [Option<Point3>], index: usize, candidate: Point3) {
    controls[index] = Some(match controls[index] {
        Some(existing) => midpoint(existing, candidate),
        None => candidate,
    });
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

fn ribbon_derivative(ribbon: &PreparedExactBoundaryRibbon, parameter: f64) -> Vec3 {
    ribbon
        .samples
        .iter()
        .min_by(|a, b| {
            (a.parameter - parameter)
                .abs()
                .partial_cmp(&(b.parameter - parameter).abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|sample| sample.cross_derivative)
        .unwrap_or(Vec3::zero())
}

fn has_gross_foldover(surface: &NurbsSurface, budget: ExactBoundarySolveBudget) -> bool {
    for u in [0.25, 0.5, 0.75] {
        for v in [0.25, 0.5, 0.75] {
            if surface.normal_at(u, v).length() < budget.min_normal_length {
                return true;
            }
        }
    }
    false
}

fn reverse_curve(curve: &NurbsCurve) -> NurbsCurve {
    let first = *curve.knots.first().unwrap_or(&0.0);
    let last = *curve.knots.last().unwrap_or(&1.0);
    let mut knots: Vec<f64> = curve
        .knots
        .iter()
        .rev()
        .map(|knot| first + last - *knot)
        .collect();
    knots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mut control_points = curve.control_points.clone();
    control_points.reverse();
    NurbsCurve {
        degree: curve.degree,
        knots,
        control_points,
        weights: curve.weights.iter().rev().copied().collect(),
    }
}

fn index(i: usize, j: usize, num_v: usize) -> usize {
    i * num_v + j
}

fn lerp(a: Point3, b: Point3, t: f64) -> Point3 {
    a + (b - a) * t
}

fn blend_point(a: Point3, b: Point3, t: f64) -> Point3 {
    lerp(a, b, t)
}

fn bilinear_point(p00: Point3, p10: Point3, p01: Point3, p11: Point3, u: f64, v: f64) -> Point3 {
    let a = blend_point(p00, p10, u);
    let b = blend_point(p01, p11, u);
    blend_point(a, b, v)
}

fn add_vec(a: Point3, b: Point3) -> Point3 {
    Point3::new(a.x + b.x, a.y + b.y, a.z + b.z)
}

fn midpoint(a: Point3, b: Point3) -> Point3 {
    a + (b - a) * 0.5
}
