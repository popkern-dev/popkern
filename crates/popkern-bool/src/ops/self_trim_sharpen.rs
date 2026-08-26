//! Intersection-vertex sharpening for the self-trim (PK-327).
//!
//! The splitter places boundary pre-split vertices at CLIP-estimated
//! parameters (UV-sampled crossings), which land within the boolean
//! era's accuracy class (~1e-3) but off the exact intersection curve
//! by more than the kernel's linear tolerance — the split body then
//! fails loop-closure validation. This is the numerical sharpening
//! step of the distance-map trimming scheme (Seong/Elber/Kim 2006
//! §5): snap each intersection-edge endpoint vertex onto the exact
//! curve, then re-trim the curve ranges of every affected edge
//! (the intersection edges themselves included — a closed circle's
//! clip-estimated span rarely ends exactly on its anchor) to the
//! final vertex positions.

use std::collections::BTreeSet;

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, VertexHandle};

use crate::error::BoolError;
use crate::split::FragmentMap;

/// Snap intersection-edge endpoints onto their exact curves, then
/// re-trim affected edges to the moved vertices. Moves are capped at
/// `10⁴ × linear` — a vertex further off its curve is a split defect,
/// not clip slop.
pub(super) fn sharpen_intersection_vertices(
    body: &mut BRep,
    fragments: &FragmentMap,
    tolerance: &Tolerance,
) -> Result<(), BoolError> {
    let tol = tolerance.linear();
    let cap = tol * 1.0e4;

    // Pass 1: snap every intersection-edge endpoint vertex onto its
    // exact curve.
    let mut touched: BTreeSet<VertexHandle> = BTreeSet::new();
    for &eh in &fragments.intersection_edges {
        let Some(edge) = body.edges.get(eh) else {
            continue;
        };
        let Some(curve) = edge.curve.as_ref() else {
            continue;
        };
        let geometry = curve.geometry.clone();
        for vh in edge_endpoints(body, eh) {
            if !touched.insert(vh) {
                continue;
            }
            let p = body.vertices[vh].position;
            let Some(t) = geometry.closest_point(&p) else {
                continue;
            };
            let q = geometry.point_at(t);
            let gap = p.distance_to(&q);
            if gap > cap {
                return Err(BoolError::SelfTrimFailed(format!(
                    "split vertex {vh:?} sits {gap:.3e} off its intersection curve \
                     (cap {cap:.3e}) — split defect, not clip slop"
                )));
            }
            if gap > tol {
                body.vertices[vh].position = q;
            }
        }
    }

    // Pass 2: re-trim the curve range of every edge touching a
    // sharpened vertex (intersection edges included) so curve-eval
    // endpoints coincide with the endpoint vertices.
    let mut edges: BTreeSet<EdgeHandle> = fragments.intersection_edges.iter().copied().collect();
    for (_, he) in body.half_edges.iter() {
        if touched.contains(&he.origin) {
            edges.insert(he.edge);
        }
    }
    for eh in edges {
        retrim_edge_to_vertices(body, eh, tol)?;
    }
    Ok(())
}

/// The (start, end) vertices of an edge's canonical half-edge.
fn edge_endpoints(body: &BRep, eh: EdgeHandle) -> Vec<VertexHandle> {
    let Some(edge) = body.edges.get(eh) else {
        return Vec::new();
    };
    let Some(he) = body.half_edges.get(edge.half_edge) else {
        return Vec::new();
    };
    let mut out = vec![he.origin];
    if let Some(twin) = body.half_edges.get(he.twin) {
        out.push(twin.origin);
    } else if let Some(next) = body.half_edges.get(he.next) {
        out.push(next.origin);
    }
    out
}

/// Clamp `eh`'s curve range so its evaluations at `t_start`/`t_end`
/// coincide with the edge's endpoint vertex positions.
fn retrim_edge_to_vertices(body: &mut BRep, eh: EdgeHandle, tol: f64) -> Result<(), BoolError> {
    let Some(edge) = body.edges.get(eh) else {
        return Ok(());
    };
    let Some(curve) = edge.curve.as_ref() else {
        return Ok(());
    };
    let geometry = curve.geometry.clone();
    let (old_start, old_end) = (curve.t_start, curve.t_end);
    let ends = edge_endpoints(body, eh);
    if ends.is_empty() {
        return Ok(());
    }
    let v_start = ends[0];
    let v_end = *ends.last().unwrap_or(&v_start);
    let p_start: Point3 = body.vertices[v_start].position;
    let p_end: Point3 = body.vertices[v_end].position;

    let project = |p: &Point3| -> Result<f64, BoolError> {
        let Some(t) = geometry.closest_point(p) else {
            return Err(BoolError::SelfTrimFailed(format!(
                "edge {eh:?}: endpoint projection failed during sharpening"
            )));
        };
        let residual = geometry.point_at(t).distance_to(p);
        if residual > tol * 100.0 {
            return Err(BoolError::SelfTrimFailed(format!(
                "edge {eh:?} cannot re-trim to its sharpened vertex \
                 (residual {residual:.3e})"
            )));
        }
        Ok(t)
    };

    if v_start == v_end || p_start.distance_to(&p_end) <= tol {
        // Closed edge: re-anchor, preserve the signed span.
        let already_ok = geometry.point_at(old_start).distance_to(&p_start) <= tol
            && geometry.point_at(old_end).distance_to(&p_end) <= tol;
        if already_ok {
            return Ok(());
        }
        let t_new = project(&p_start)?;
        let span = old_end - old_start;
        if let Some(c) = body.edges[eh].curve.as_mut() {
            c.t_start = t_new;
            c.t_end = t_new + span;
        }
        return Ok(());
    }

    // Open edge: match ends to vertices by proximity of the CURRENT
    // evaluations (the range orientation is already established).
    let eval_start = geometry.point_at(old_start);
    let eval_end = geometry.point_at(old_end);
    let (target_start, target_end) = if eval_start.distance_to(&p_start)
        + eval_end.distance_to(&p_end)
        <= eval_start.distance_to(&p_end) + eval_end.distance_to(&p_start)
    {
        (p_start, p_end)
    } else {
        (p_end, p_start)
    };
    if eval_start.distance_to(&target_start) > tol {
        let t = project(&target_start)?;
        if let Some(c) = body.edges[eh].curve.as_mut() {
            c.t_start = t;
        }
    }
    if eval_end.distance_to(&target_end) > tol {
        let t = project(&target_end)?;
        if let Some(c) = body.edges[eh].curve.as_mut() {
            c.t_end = t;
        }
    }
    Ok(())
}
