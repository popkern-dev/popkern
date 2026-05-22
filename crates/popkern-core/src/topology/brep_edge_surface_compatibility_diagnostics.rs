//! Structured diagnostics for edge/surface compatibility repair.

use crate::geometry::edge_triple::{BoundarySide, EdgeHint, SurfaceSide};
use crate::geometry::surface::SurfaceGeometry;
use crate::topology::{SharedTrimEdgeFailureReason, SurfaceBoundaryHintDeclineReason};

/// Surface class observed by edge/surface compatibility repair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BrepEdgeSurfaceKind {
    /// Analytic plane.
    Plane,
    /// Analytic cylinder.
    Cylinder,
    /// Analytic cone.
    Cone,
    /// Analytic sphere.
    Sphere,
    /// Analytic torus.
    Torus,
    /// NURBS surface.
    Nurbs,
    /// Gregory surface.
    Gregory,
}

/// Boundary hint considered during edge/surface compatibility repair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BrepEdgeSurfaceBoundaryHintEvidence {
    /// Which side owned the inferred boundary.
    pub side: SurfaceSide,
    /// Boundary selected on the owning surface.
    pub boundary: BoundarySide,
}

/// Authored boundary evidence that was rejected before repair trusted it.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepEdgeSurfaceBoundaryDeclineEvidence {
    /// Which adjacent side owned the rejected hint.
    pub side: SurfaceSide,
    /// Machine-readable decline reason.
    pub reason: SurfaceBoundaryHintDeclineReason,
    /// Best alternative boundary measured from the pcurve, when available.
    pub best_boundary: Option<BoundarySide>,
    /// Best sampled deviation from any canonical boundary.
    pub best_boundary_deviation: f64,
    /// Maximum topological endpoint gap to the claimed boundary trace.
    pub endpoint_gap: f64,
}

/// Per-edge reconstruction evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepEdgeSurfaceRepairEvidence {
    /// Surface class on the first adjacent face use.
    pub surface_a: Option<BrepEdgeSurfaceKind>,
    /// Surface class on the second adjacent face use.
    pub surface_b: Option<BrepEdgeSurfaceKind>,
    /// Boundary hint inferred or authored for this edge.
    pub boundary_hint: Option<BrepEdgeSurfaceBoundaryHintEvidence>,
    /// Authored boundary hints rejected before repair trusted them.
    pub boundary_declines: Vec<BrepEdgeSurfaceBoundaryDeclineEvidence>,
    /// Whether boundary-owned reconstruction was attempted.
    pub attempted_boundary_reconstruction: bool,
    /// Whether generic surface-intersection reconstruction was attempted.
    pub attempted_surface_intersection: bool,
    /// Whether unanchored intersection plus vertex re-anchoring was attempted.
    pub attempted_anchor_repair: bool,
    /// Maximum vertex motion proposed or accepted by anchor repair.
    pub anchor_motion: f64,
    /// Maximum vertex motion allowed by the caller's local repair policy.
    pub max_allowed_anchor_motion: f64,
    /// Maximum endpoint gap this anchor move would introduce on incident edges.
    pub max_incident_endpoint_gap: f64,
    /// Diagnostic label for the boundary attempt failure, if any.
    pub boundary_failure: Option<String>,
    /// Diagnostic label for the intersection attempt failure, if any.
    pub intersection_failure: Option<String>,
    /// Shared trim-edge validation failure from the atomic commit primitive, if any.
    pub shared_trim_failure: Option<SharedTrimEdgeFailureReason>,
    /// Shared trim-edge candidate curve/surface deviation at commit validation.
    pub shared_trim_curve_surface_deviation: f64,
    /// Shared trim-edge candidate side-to-side deviation at commit validation.
    pub shared_trim_side_agreement_deviation: f64,
    /// Shared trim-edge candidate endpoint gap at commit validation.
    pub shared_trim_endpoint_gap: f64,
    /// Candidate curve residual against adjacent side A.
    pub residual_deviation_a: f64,
    /// Candidate curve residual against adjacent side B.
    pub residual_deviation_b: f64,
    /// Stable diagnostic codes from the scoped repair knowledge graph.
    pub repair_knowledge_diagnostics: Vec<String>,
}

impl BrepEdgeSurfaceRepairEvidence {
    /// Evidence for declines before adjacent surfaces are known.
    pub fn unavailable() -> Self {
        Self {
            surface_a: None,
            surface_b: None,
            boundary_hint: None,
            boundary_declines: Vec::new(),
            attempted_boundary_reconstruction: false,
            attempted_surface_intersection: false,
            attempted_anchor_repair: false,
            anchor_motion: 0.0,
            max_allowed_anchor_motion: 0.0,
            max_incident_endpoint_gap: 0.0,
            boundary_failure: None,
            intersection_failure: None,
            shared_trim_failure: None,
            shared_trim_curve_surface_deviation: 0.0,
            shared_trim_side_agreement_deviation: 0.0,
            shared_trim_endpoint_gap: 0.0,
            residual_deviation_a: 0.0,
            residual_deviation_b: 0.0,
            repair_knowledge_diagnostics: Vec::new(),
        }
    }

    /// Evidence initialized from known adjacent surfaces and a possible hint.
    pub fn from_surfaces(
        surface_a: &SurfaceGeometry,
        surface_b: &SurfaceGeometry,
        hint: Option<&EdgeHint>,
    ) -> Self {
        Self {
            surface_a: Some(surface_kind(surface_a)),
            surface_b: Some(surface_kind(surface_b)),
            boundary_hint: hint.and_then(boundary_hint_evidence),
            boundary_declines: Vec::new(),
            attempted_boundary_reconstruction: false,
            attempted_surface_intersection: false,
            attempted_anchor_repair: false,
            anchor_motion: 0.0,
            max_allowed_anchor_motion: 0.0,
            max_incident_endpoint_gap: 0.0,
            boundary_failure: None,
            intersection_failure: None,
            shared_trim_failure: None,
            shared_trim_curve_surface_deviation: 0.0,
            shared_trim_side_agreement_deviation: 0.0,
            shared_trim_endpoint_gap: 0.0,
            residual_deviation_a: 0.0,
            residual_deviation_b: 0.0,
            repair_knowledge_diagnostics: Vec::new(),
        }
    }

    /// Evidence initialized for a one-sided surfaced trim whose twin face has
    /// no repairable surface. This is still a valid repair class for sheet or
    /// residual construction topology: the repaired edge must lie on the one
    /// surfaced owner and preserve the existing topological endpoints.
    pub fn from_single_surface(surface: &SurfaceGeometry) -> Self {
        Self {
            surface_a: Some(surface_kind(surface)),
            surface_b: None,
            boundary_hint: None,
            boundary_declines: Vec::new(),
            attempted_boundary_reconstruction: true,
            attempted_surface_intersection: false,
            attempted_anchor_repair: false,
            anchor_motion: 0.0,
            max_allowed_anchor_motion: 0.0,
            max_incident_endpoint_gap: 0.0,
            boundary_failure: None,
            intersection_failure: None,
            shared_trim_failure: None,
            shared_trim_curve_surface_deviation: 0.0,
            shared_trim_side_agreement_deviation: 0.0,
            shared_trim_endpoint_gap: 0.0,
            residual_deviation_a: 0.0,
            residual_deviation_b: 0.0,
            repair_knowledge_diagnostics: Vec::new(),
        }
    }
}

/// Classify a surface for diagnostics without exposing implementation data.
pub fn surface_kind(surface: &SurfaceGeometry) -> BrepEdgeSurfaceKind {
    match surface {
        SurfaceGeometry::Plane { .. } => BrepEdgeSurfaceKind::Plane,
        SurfaceGeometry::Cylinder { .. } => BrepEdgeSurfaceKind::Cylinder,
        SurfaceGeometry::Cone { .. } => BrepEdgeSurfaceKind::Cone,
        SurfaceGeometry::Sphere { .. } => BrepEdgeSurfaceKind::Sphere,
        SurfaceGeometry::Torus { .. } => BrepEdgeSurfaceKind::Torus,
        SurfaceGeometry::Nurbs(_) => BrepEdgeSurfaceKind::Nurbs,
        SurfaceGeometry::Gregory(_) => BrepEdgeSurfaceKind::Gregory,
    }
}

fn boundary_hint_evidence(hint: &EdgeHint) -> Option<BrepEdgeSurfaceBoundaryHintEvidence> {
    match hint {
        EdgeHint::SurfaceBoundary { side, boundary } => Some(BrepEdgeSurfaceBoundaryHintEvidence {
            side: *side,
            boundary: *boundary,
        }),
        _ => None,
    }
}
