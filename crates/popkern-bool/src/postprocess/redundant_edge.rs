//! Redundant edge detection and removal.
//!
//! After a boolean, some edges separate two faces that lie on the same
//! surface with the same orientation. These edges are geometrically
//! unnecessary and can be removed to simplify the result.

use popkern_core::geometry::intersection::coincident_analytic::are_surfaces_coincident_analytic;
use popkern_core::geometry::surface::SurfaceOnFace;
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
    find_redundant_edges_in(brep, tolerance, &crate::serve_region::ServeRegion::Whole)
}

/// [`find_redundant_edges`] over the edges `region` admits (both faces
/// in it).
pub(crate) fn find_redundant_edges_in(
    brep: &BRep,
    tolerance: &Tolerance,
    region: &crate::serve_region::ServeRegion,
) -> RedundantEdgeReport {
    let mut report = RedundantEdgeReport::default();

    for (eh, _) in &brep.edges {
        if !region.admits_edge(brep, eh) {
            continue;
        }
        report.edges_analyzed += 1;
        if edge_redundant_now(brep, eh, tolerance) {
            report.redundant_edges.push(eh);
        }
    }

    report
}

/// The full redundancy predicate against the CURRENT topology. The
/// KEF loop re-evaluates this immediately before each dissolution:
/// the report is computed once, and earlier KEFs change which edges
/// two faces share — a stale single-share verdict dissolved one edge
/// of a (by then) multi-edge shared boundary, and the survivor came
/// out SELF-TWINNED inside the merged loop (measured on the torus ×
/// cylinder union: 20 self-twinned grid edges, kemr −1 half-edge
/// delta, STEP demoted to OPEN_SHELL).
pub(crate) fn edge_redundant_now(brep: &BRep, eh: EdgeHandle, tolerance: &Tolerance) -> bool {
    (|| -> Option<bool> {
        let edge = brep.edges.get(eh)?;

        let he = brep.half_edges.get(edge.half_edge)?;
        let twin = brep.half_edges.get(he.twin)?;
        // A self-twinned half-edge is already residue, never a
        // dissolution candidate.
        if he.twin == edge.half_edge {
            return Some(false);
        }
        let loop_a = brep.loops.get(he.loop_)?;
        let loop_b = brep.loops.get(twin.loop_)?;

        let face_a = loop_a.face;
        let face_b = loop_b.face;

        // Skip seam edges (face_a == face_b)
        if face_a == face_b {
            return Some(false);
        }

        let surf_a = brep.faces.get(face_a).and_then(|f| f.surface.as_ref())?;
        let surf_b = brep.faces.get(face_b).and_then(|f| f.surface.as_ref())?;

        // Check same surface (analytic comparison only; the locus,
        // whatever its parameterisation or normal sign)
        if !are_surfaces_coincident_analytic(&surf_a.geometry, &surf_b.geometry, tolerance) {
            return Some(false);
        }

        if !outward_normals_agree(brep, eh, surf_a, surf_b)? {
            return Some(false);
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
        Some(shared_count == 1)
    })()
    .unwrap_or(false)
}

/// Same ORIENTED surface: the two faces' outward normals agree at edge
/// `eh`'s midpoint. The stored flags alone are not a comparison: a cap
/// with normal +z marked reversed IS a face with normal −z marked
/// forward (measured, PK-750: a rotated rod's bottom cap stayed split
/// from a flush box's coplanar bottom along their chord while the top
/// merged; PK-778: the same flag test in the shared-CHAIN path kept the
/// caps split wherever the rod's seam vertex divided the shared rim).
pub(crate) fn outward_normals_agree(
    brep: &BRep,
    eh: EdgeHandle,
    surf_a: &SurfaceOnFace,
    surf_b: &SurfaceOnFace,
) -> Option<bool> {
    use popkern_core::geometry::curve::CurveEval as _;
    use popkern_core::geometry::surface::SurfaceEval as _;
    let c = brep.edges.get(eh)?.curve.as_ref()?;
    let mid = c.geometry.point_at(0.5 * (c.t_start + c.t_end));
    let outward = |s: &SurfaceOnFace| {
        let (u, v) = s.geometry.closest_point(&mid)?;
        let n = s.geometry.normal_at(u, v);
        Some(if s.orientation { n } else { n * -1.0 })
    };
    Some(outward(surf_a)?.dot(&outward(surf_b)?) > 0.0)
}
