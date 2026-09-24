//! Corner-resolution postcondition (PK-516).
//!
//! A multi-fillet corner that the vertex phase leaves WITHOUT a
//! corner blend — `VertexMode::Mitred` requested, the smooth vertex
//! blend falling back, or an N >= 5 all-or-nothing rollback — is
//! served by whatever the face-trim pipeline left at the corner.
//! That residue was laid out for the corner blend to replace: the
//! strips run on past their true section toward the original apex,
//! bounded by straight scaffold edges that lie on the support planes
//! but NOT on the strip surfaces (measured on the pentagonal pyramid,
//! r = 0.16: ten rim edges 1.87e-2 off their NURBS strips). The body
//! is structurally sound, so the fast `validate()` (which excludes
//! the check-9 curve-surface sweep by design) passes, and before this
//! gate `fillet_edges` returned it as `Ok` — an invalid body that only
//! the part layer's `validate_all` caught, far from the cause.
//!
//! Every unresolved corner is therefore recorded with the face star
//! it owned when the corner machinery gave up, and the pipeline exit
//! measures curve-surface consistency (check 9) over the edges of
//! those faces AFTER every tail pass has run — a corner that a later
//! pass closes exactly (the blend-seam corner primitive at valence-2
//! corners) serves, one whose residue is still off-surface declines
//! the whole fillet transactionally with `VertexBlendFailed`, naming
//! the vertex, why its blend was not built, and the worst edge. The
//! gate measures only the corner territory: it neither pays for nor
//! fails on a whole-body sweep of geometry the corner did not own.
//! Tolerance is each edge's own band (check 9's contract), not a
//! gate-local constant.

use popkern_core::geometry::Point3;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

/// Why a corner has no corner blend.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum UnresolvedCornerCause {
    /// The caller asked for `VertexMode::Mitred`.
    MitredRequested,
    /// The smooth vertex blend declined; the reason it gave.
    SmoothBlendDeclined(String),
    /// An N >= 5 corner rolled back because no route produced an
    /// accepted N-patch tiling; each route's final decline verdict
    /// (`stage: detail`), in dispatch order.
    AllOrNothingRollback(Vec<String>),
}

impl std::fmt::Display for UnresolvedCornerCause {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::MitredRequested => write!(f, "mitred corner requested"),
            Self::SmoothBlendDeclined(reason) => {
                write!(f, "smooth vertex blend declined: {reason}")
            }
            Self::AllOrNothingRollback(route_declines) => {
                write!(
                    f,
                    "no N-patch corner tiling accepted (all-or-nothing rollback)"
                )?;
                if !route_declines.is_empty() {
                    write!(f, " — {}", route_declines.join("; "))?;
                }
                Ok(())
            }
        }
    }
}

/// A corner served without a corner blend, with the territory the
/// exit gate measures.
#[derive(Debug, Clone)]
pub struct UnresolvedCorner {
    /// The corner vertex (may no longer exist in the served body).
    pub vertex: VertexHandle,
    /// Its position when the corner machinery gave up.
    pub position: Point3,
    /// A filleted edge meeting the corner, for failure attribution.
    pub edge: EdgeHandle,
    /// Faces incident to the vertex at that moment.
    pub faces: Vec<FaceHandle>,
    /// Why no corner blend was built.
    pub cause: UnresolvedCornerCause,
}

/// Record an unresolved corner: snapshot the face star of `vertex`
/// in `brep`. Returns `None` when the vertex is already gone (the
/// corner machinery consumed it — nothing unresolved remains there).
pub(crate) fn record_unresolved_corner(
    brep: &BRep,
    vertex: VertexHandle,
    edge: EdgeHandle,
    cause: UnresolvedCornerCause,
) -> Option<UnresolvedCorner> {
    let position = brep.vertices.get(vertex)?.position;
    let mut faces: Vec<FaceHandle> = Vec::new();
    for he in brep.half_edges.values() {
        if he.origin != vertex {
            continue;
        }
        let Some(face) = brep.loops.get(he.loop_).map(|l| l.face) else {
            continue;
        };
        if !faces.contains(&face) {
            faces.push(face);
        }
    }
    Some(UnresolvedCorner {
        vertex,
        position,
        edge,
        faces,
        cause,
    })
}

/// The final decline verdict of each corner cap route, from the
/// structured apex diagnostics: acceptance-gate declines and
/// construction errors (`stage: detail`). Intermediate candidate
/// facts are left to the repair report.
pub(crate) fn route_decline_verdicts(
    diagnostics: &[crate::types::ApexConstructionDiagnostic],
) -> Vec<String> {
    diagnostics
        .iter()
        .filter_map(|diagnostic| match diagnostic {
            crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(d)
                if d.code.ends_with(".declined")
                    && (d.stage.ends_with("cap_support_acceptance")
                        || d.stage.ends_with("construction_error")
                        || d.stage.ends_with("self_intersection_overlap")) =>
            {
                Some(format!(
                    "{}: {}",
                    d.stage,
                    d.detail_code.as_deref().unwrap_or("declined")
                ))
            }
            _ => None,
        })
        .collect()
}

/// Every multi-fillet corner of a `VertexMode::Mitred` request: the
/// vertex phase builds no corner blend at any of them.
pub(crate) fn requested_mitre_corners(
    brep: &BRep,
    edge_geometries: &[super::EdgeGeometry],
) -> Vec<UnresolvedCorner> {
    crate::analysis::vertex_incidences(
        edge_geometries
            .iter()
            .map(|eg| (eg.analysis.vertex_start, eg.analysis.vertex_end)),
    )
    .into_iter()
    .filter(|(_, edge_indices)| edge_indices.len() >= 2)
    .filter_map(|(vertex, edge_indices)| {
        record_unresolved_corner(
            brep,
            vertex,
            edge_geometries[edge_indices[0]].analysis.edge,
            UnresolvedCornerCause::MitredRequested,
        )
    })
    .collect()
}

/// Edges of every surviving face in `faces` (outer and inner loops),
/// deduplicated, in first-seen order.
fn territory_edges(brep: &BRep, faces: &[FaceHandle]) -> Vec<EdgeHandle> {
    let mut edges: Vec<EdgeHandle> = Vec::new();
    for &face in faces {
        let Some(face_data) = brep.faces.get(face) else {
            continue;
        };
        for &loop_handle in std::iter::once(&face_data.outer_loop).chain(&face_data.inner_loops) {
            let Some(start) = brep.loops.get(loop_handle).map(|l| l.half_edge) else {
                continue;
            };
            let mut cursor = start;
            for _ in 0..=brep.half_edges.len() {
                let Some(he) = brep.half_edges.get(cursor) else {
                    break;
                };
                if !edges.contains(&he.edge) {
                    edges.push(he.edge);
                }
                cursor = he.next;
                if cursor == start {
                    break;
                }
            }
        }
    }
    edges
}

/// The measured violation that declines an unresolved corner.
#[derive(Debug, Clone)]
pub(crate) struct CornerResolutionViolation {
    /// The corner.
    pub corner: UnresolvedCorner,
    /// Its worst off-surface edge.
    pub worst_edge: EdgeHandle,
    /// That edge's worst curve-surface deviation.
    pub worst_deviation: f64,
    /// How many territory edges fail check 9.
    pub failing_edges: usize,
}

impl CornerResolutionViolation {
    /// Human-readable decline reason for `VertexBlendFailed`.
    pub(crate) fn reason(&self) -> String {
        format!(
            "corner left without a corner blend ({}) and its territory fails curve-surface \
             consistency: {} edge(s) off their faces, worst {:?} at {:.3e}",
            self.corner.cause, self.failing_edges, self.worst_edge, self.worst_deviation
        )
    }
}

/// Measure every unresolved corner's territory on the served body;
/// the first corner (recording order) with any check-9 violation is
/// returned.
pub(crate) fn first_corner_resolution_violation(
    brep: &BRep,
    corners: &[UnresolvedCorner],
) -> Option<CornerResolutionViolation> {
    for corner in corners {
        let errors = brep.validate_curve_surface_edges(&territory_edges(brep, &corner.faces));
        let worst = errors
            .iter()
            .filter_map(|error| match error {
                popkern_core::error::ValidationError::CurveSurfaceDeviation {
                    edge,
                    max_deviation,
                } => Some((*edge, *max_deviation)),
                _ => None,
            })
            .max_by(|a, b| a.1.total_cmp(&b.1));
        if let Some((worst_edge, worst_deviation)) = worst {
            return Some(CornerResolutionViolation {
                corner: corner.clone(),
                worst_edge,
                worst_deviation,
                failing_edges: errors.len(),
            });
        }
    }
    None
}
