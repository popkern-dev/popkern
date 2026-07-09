//! Redundant edge detection and removal.
//!
//! After a boolean, some edges separate two faces that lie on the same
//! surface with the same orientation. These edges are geometrically
//! unnecessary and can be removed to simplify the result.

use popkern_core::geometry::intersection::coincident_analytic::are_surfaces_coincident_analytic;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle};

/// Result of redundant edge analysis.
#[derive(Debug, Clone, Default)]
pub struct RedundantEdgeReport {
    /// Edge handles that are redundant (can be removed).
    pub redundant_edges: Vec<EdgeHandle>,
    /// Total edges analyzed.
    pub edges_analyzed: usize,
}

/// Find edges between co-surfacial, same-orientation faces.
///
/// An edge is redundant if:
/// 1. Its two incident faces have the same underlying surface
/// 2. The faces have the same orientation
/// 3. Both faces lie on an analytic surface (Plane, Cylinder, etc.)
///
/// NURBS faces are excluded from simplification because co-surfacial
/// detection for NURBS is approximate and unreliable.
pub fn find_redundant_edges(brep: &BRep, tolerance: &Tolerance) -> RedundantEdgeReport {
    let mut report = RedundantEdgeReport::default();

    for (eh, edge) in &brep.edges {
        report.edges_analyzed += 1;

        let he = match brep.half_edges.get(edge.half_edge) {
            Some(h) => h,
            None => continue,
        };
        let twin = match brep.half_edges.get(he.twin) {
            Some(t) => t,
            None => continue,
        };

        let loop_a = match brep.loops.get(he.loop_) {
            Some(l) => l,
            None => continue,
        };
        let loop_b = match brep.loops.get(twin.loop_) {
            Some(l) => l,
            None => continue,
        };

        let face_a = loop_a.face;
        let face_b = loop_b.face;

        // Skip seam edges (face_a == face_b)
        if face_a == face_b {
            continue;
        }

        let surf_a = match brep.faces.get(face_a).and_then(|f| f.surface.as_ref()) {
            Some(s) => s,
            None => continue,
        };
        let surf_b = match brep.faces.get(face_b).and_then(|f| f.surface.as_ref()) {
            Some(s) => s,
            None => continue,
        };

        // Check same orientation
        if surf_a.orientation != surf_b.orientation {
            continue;
        }

        // Check same surface (analytic comparison only)
        if !are_surfaces_coincident_analytic(&surf_a.geometry, &surf_b.geometry, tolerance) {
            continue;
        }

        // KEF on this edge is only sound when it is the ONLY edge the
        // two faces share: dissolving one edge of a multi-edge shared
        // boundary leaves the survivors traversed twice by the merged
        // loop — a zero-width slit that passes Euler and manifold
        // validation but cannot be triangulated (measured on a flush
        // coplanar union: an arm wall and the gusset cap shared both
        // the slope record and the corner column; the merged face
        // lost 70% of the swept volume to mesh failures).
        let shared_count = brep
            .edges
            .iter()
            .filter(|(_, other)| {
                let Some(ohe) = brep.half_edges.get(other.half_edge) else {
                    return false;
                };
                let Some(otwin) = brep.half_edges.get(ohe.twin) else {
                    return false;
                };
                let (Some(ola), Some(olb)) =
                    (brep.loops.get(ohe.loop_), brep.loops.get(otwin.loop_))
                else {
                    return false;
                };
                (ola.face == face_a && olb.face == face_b)
                    || (ola.face == face_b && olb.face == face_a)
            })
            .count();
        if shared_count == 1 {
            report.redundant_edges.push(eh);
        }
    }

    report
}
