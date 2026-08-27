//! Post-corner seam re-authoring (PK-344). Split from
//! `edge_recompute.rs` (500-line ratchet).

use crate::error::ShellError;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle};

use super::edge_curve_dispatch::compute_edge_curve;

/// Re-author a SAME-SURFACE (seam) edge whose curve endpoints no
/// longer match its vertices (PK-344): gap-corner anchors carry
/// interim positions during edge recompute (the corner phase places
/// them later), so a seam ending at such an anchor declines its
/// recompute and stays stale — the corner construction calls this
/// AFTER moving the vertices to their true targets. Gated on measured
/// endpoint disagreement so already-consistent seams are untouched
/// (bit-stability of serving bodies).
///
/// # Errors
///
/// Propagates topology lookups; a declined ladder leaves the curve
/// unchanged (the serve gate owns the verdict).
pub(crate) fn reauthor_stale_same_surface_edge(
    inner: &mut BRep,
    edge_handle: EdgeHandle,
    tolerance: &Tolerance,
) -> Result<(), ShellError> {
    let Some(edge) = inner.edges.get(edge_handle) else {
        return Ok(());
    };
    let Some(he) = inner.half_edges.get(edge.half_edge) else {
        return Ok(());
    };
    let Some(twin) = inner.half_edges.get(he.twin) else {
        return Ok(());
    };
    let (l1, l2) = (he.loop_, twin.loop_);
    let (Some(f1), Some(f2)) = (
        inner.loops.get(l1).map(|l| l.face),
        inner.loops.get(l2).map(|l| l.face),
    ) else {
        return Ok(());
    };
    if f1 != f2 {
        return Ok(()); // not a seam
    }
    let Some(surface) = inner
        .faces
        .get(f1)
        .and_then(|f| f.surface.as_ref())
        .map(|s| s.geometry.clone())
    else {
        return Ok(());
    };
    let (Some(v0), Some(v1)) = (
        inner.vertices.get(he.origin).map(|v| v.position),
        inner.vertices.get(twin.origin).map(|v| v.position),
    ) else {
        return Ok(());
    };
    let stale = match edge.curve.as_ref() {
        Some(c) => {
            use popkern_core::geometry::curve::CurveEval;
            let a = c.geometry.point_at(c.t_start);
            let b = c.geometry.point_at(c.t_end);
            let band = tolerance.linear() * 10.0;
            !((a.distance_to(&v0) < band && b.distance_to(&v1) < band)
                || (a.distance_to(&v1) < band && b.distance_to(&v0) < band))
        }
        None => false,
    };
    if !stale {
        return Ok(());
    }
    if let Some(curve_on_edge) =
        compute_edge_curve(inner, edge_handle, &surface, &surface, tolerance)?
    {
        inner.edges[edge_handle].curve = Some(curve_on_edge);
    }
    Ok(())
}
