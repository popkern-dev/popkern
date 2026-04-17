//! §5.3 vertex blend classifier.
//!
//! Given a vertex and the list of filleted edges incident to it,
//! [`classify_vertex_blend`] returns one of five
//! [`VertexBlendClassification`] variants matching the spec's
//! "same-side, opposite-side, mixed, tangent, free boundary"
//! taxonomy. The dispatcher uses this to pick the right corner
//! construction (rolling-ball sphere, Coons patch, Gregory patch,
//! plate minimization, or a graceful-degrade fallback).
//!
//! Public — UI layers call this before a fillet op to preview
//! which corner method will fire and warn the user about
//! potentially-hard configurations.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, HalfEdgeHandle, VertexHandle};

use crate::analysis::{Convexity, analyze_edges};

/// §5.3 corner configuration class. Drives the dispatcher's
/// choice of construction method.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VertexBlendClassification {
    /// All incident filleted edges are **convex** and the rolling
    /// ball sits on the same side of every support face. The
    /// classic cube-corner case — cleanly handled by a rolling-
    /// ball sphere for equal-radius trihedral, or by Gregory G1/G2
    /// when radii differ or N ≥ 4.
    SameSideConvex,
    /// All incident filleted edges are **concave**. Analogous to
    /// `SameSideConvex` but the rolling ball sits on the opposite
    /// side — the "pocket corner" case. Same dispatch targets.
    SameSideConcave,
    /// Some incident edges are convex, others concave. The rolling
    /// ball cannot sit uniformly; construction strategies specific
    /// to saddle-type vertices apply. This is the hardest §5.3
    /// case outside of tangent configurations.
    OppositeSide,
    /// At least one incident edge is mixed-convexity (its own
    /// classification is [`Convexity::Mixed`]), or the convexity
    /// set doesn't cleanly split into convex/concave sub-groups.
    /// A UI should flag this as potentially problematic.
    Mixed,
    /// At least one incident filleted edge is **tangent** — the
    /// dihedral angle at the edge is within `tolerance.angular()`
    /// of π. The rolling-ball contact solution is undefined here;
    /// the §8.4 natural run-out or a tangent-regularization
    /// strategy applies.
    Tangent,
    /// The vertex sits on an open-shell boundary — at least one
    /// incident edge is a boundary edge (null twin) or the vertex
    /// is missing a face in its 1-ring. §8.3 cap insertion applies
    /// to the boundary side; the corner patch spans only the
    /// filleted sub-region.
    FreeBoundary,
}

/// Full classification report — includes the summary
/// [`VertexBlendClassification`] plus the raw per-edge counts a
/// UI can surface in a diagnostic panel.
#[derive(Debug, Clone)]
pub struct VertexBlendClassificationReport {
    /// The summary classification driving dispatch.
    pub classification: VertexBlendClassification,
    /// Number of filleted edges incident to the vertex.
    pub n_incident_filleted_edges: usize,
    /// Number of convex-classified filleted edges.
    pub n_convex: usize,
    /// Number of concave-classified filleted edges.
    pub n_concave: usize,
    /// Number of mixed-convexity filleted edges.
    pub n_mixed: usize,
    /// Number of tangent (dihedral ≈ π) filleted edges.
    pub n_tangent: usize,
    /// Whether any incident edge sits on an open-shell boundary
    /// (null twin) — triggers the `FreeBoundary` classification
    /// regardless of other counts.
    pub has_boundary_edge: bool,
}

/// Classify the §5.3 corner configuration at `vertex`, given the
/// filleted edges incident to it.
///
/// `filleted_edges` should be the set of edges the caller plans to
/// fillet that touch `vertex`. Edges not incident to `vertex` are
/// silently skipped; missing edges are skipped too. `tolerance`
/// drives the tangent-angle test.
///
/// Returns a [`VertexBlendClassificationReport`]. If no filleted
/// edges are incident to the vertex, returns a report with
/// `SameSideConvex` as the degenerate default and zero counts —
/// callers checking `n_incident_filleted_edges == 0` can treat
/// this as "no corner work needed."
pub fn classify_vertex_blend(
    brep: &BRep,
    vertex: VertexHandle,
    filleted_edges: &[EdgeHandle],
    tolerance: &Tolerance,
) -> VertexBlendClassificationReport {
    // Filter to edges incident to the vertex.
    let incident: Vec<EdgeHandle> = filleted_edges
        .iter()
        .copied()
        .filter(|&e| edge_is_incident_to(brep, e, vertex))
        .collect();

    if incident.is_empty() {
        return VertexBlendClassificationReport {
            classification: VertexBlendClassification::SameSideConvex,
            n_incident_filleted_edges: 0,
            n_convex: 0,
            n_concave: 0,
            n_mixed: 0,
            n_tangent: 0,
            has_boundary_edge: false,
        };
    }

    let has_boundary_edge = incident.iter().any(|&e| is_boundary_edge(brep, e));

    // Analyse each incident edge for convexity. If analysis fails
    // (e.g. seam edge, missing geometry) we count it as Mixed —
    // the dispatcher will drop to a conservative strategy.
    let analyses = analyze_edges(brep, &incident, false);

    let (mut n_convex, mut n_concave, mut n_mixed) = (0, 0, 0);
    let mut n_tangent = 0;
    match &analyses {
        Ok(analyses) => {
            for a in analyses {
                if is_edge_tangent(brep, a.edge, tolerance) {
                    n_tangent += 1;
                    continue;
                }
                match a.convexity {
                    Convexity::Convex => n_convex += 1,
                    Convexity::Concave => n_concave += 1,
                    Convexity::Mixed => n_mixed += 1,
                }
            }
        }
        Err(_) => {
            // All-or-nothing analysis failure — conservative default.
            n_mixed = incident.len();
        }
    }

    let classification = if has_boundary_edge {
        VertexBlendClassification::FreeBoundary
    } else if n_tangent > 0 {
        VertexBlendClassification::Tangent
    } else if n_mixed > 0 {
        VertexBlendClassification::Mixed
    } else if n_convex > 0 && n_concave > 0 {
        VertexBlendClassification::OppositeSide
    } else if n_convex > 0 {
        VertexBlendClassification::SameSideConvex
    } else if n_concave > 0 {
        VertexBlendClassification::SameSideConcave
    } else {
        // All indeterminate — classify as Mixed to be safe.
        VertexBlendClassification::Mixed
    };

    VertexBlendClassificationReport {
        classification,
        n_incident_filleted_edges: incident.len(),
        n_convex,
        n_concave,
        n_mixed,
        n_tangent,
        has_boundary_edge,
    }
}

/// Is `edge` incident to `vertex`?
fn edge_is_incident_to(brep: &BRep, edge: EdgeHandle, vertex: VertexHandle) -> bool {
    let Some(edata) = brep.edges.get(edge) else {
        return false;
    };
    let he_h = edata.half_edge;
    let Some(he) = brep.half_edges.get(he_h) else {
        return false;
    };
    if he.origin == vertex {
        return true;
    }
    // Check twin (manifold edges) or loop-walk next (boundary).
    let twin = he.twin;
    if twin != HalfEdgeHandle::default()
        && let Some(t) = brep.half_edges.get(twin)
        && t.origin == vertex
    {
        return true;
    }
    // Boundary edge — derive second endpoint from he.next.
    if let Some(n) = brep.half_edges.get(he.next)
        && n.origin == vertex
    {
        return true;
    }
    false
}

/// Is `edge` a boundary edge (null twin)?
fn is_boundary_edge(brep: &BRep, edge: EdgeHandle) -> bool {
    brep.boundary_edges.contains(&edge)
        || brep
            .edges
            .get(edge)
            .and_then(|e| brep.half_edges.get(e.half_edge))
            .is_some_and(|he| he.twin == HalfEdgeHandle::default())
}

/// Is the dihedral angle at `edge` within `tolerance.angular()` of
/// π (tangent faces)? Uses the surface normals at the edge midpoint.
/// Returns `false` for edges without a curve or for seam edges.
fn is_edge_tangent(brep: &BRep, edge: EdgeHandle, tolerance: &Tolerance) -> bool {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::SurfaceEval;

    let Some(edata) = brep.edges.get(edge) else {
        return false;
    };
    let Some(curve) = edata.curve.as_ref() else {
        return false;
    };
    let he_h = edata.half_edge;
    let Some(he) = brep.half_edges.get(he_h) else {
        return false;
    };
    let twin = he.twin;
    if twin == HalfEdgeHandle::default() {
        return false; // boundary — no two-sided dihedral
    }
    let Some(t_he) = brep.half_edges.get(twin) else {
        return false;
    };
    let face_a = brep.loops[he.loop_].face;
    let face_b = brep.loops[t_he.loop_].face;
    if face_a == face_b {
        return false; // seam edge
    }
    let Some(sa) = brep.faces.get(face_a).and_then(|f| f.surface.as_ref()) else {
        return false;
    };
    let Some(sb) = brep.faces.get(face_b).and_then(|f| f.surface.as_ref()) else {
        return false;
    };
    let t_mid = 0.5 * (curve.t_start + curve.t_end);
    let p = curve.geometry.point_at(t_mid);
    let Some((u_a, v_a)) = sa.geometry.closest_point(&p) else {
        return false;
    };
    let Some((u_b, v_b)) = sb.geometry.closest_point(&p) else {
        return false;
    };
    let n_a = sa.geometry.normal_at(u_a, v_a).normalized();
    let n_b = sb.geometry.normal_at(u_b, v_b).normalized();
    let (Some(n_a), Some(n_b)) = (n_a, n_b) else {
        return false;
    };
    // Dihedral at π means n_a · n_b ≈ -1 (normals opposite because
    // faces fold back to near-coplanar). The angular tolerance in
    // cosine terms: cos(π - tol) = -cos(tol) ≈ -(1 - tol²/2).
    let cos_tol = tolerance.angular().cos();
    n_a.dot(&n_b) < -cos_tol
}

#[cfg(test)]
#[path = "vertex_classifier_tests.rs"]
mod tests;
