//! the anchor-split stage anchor insertion: vertex reuse, straight and curved edge
//! splitting, and the curved-support interior spur (PK-141). Split
//! from `helpers.rs` (line ratchet).

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

use popkern_core::geometry::curve::CurveEval;

/// Insert a new vertex at `position` in the mutable BRep by finding
/// the descendant sub-edge currently containing it and splitting it.
/// Appends `(old_edge, edge_a, edge_b)` to `edge_splits` so the
/// caller can maintain a fillet-edge → descendants map.
///
/// If a vertex already exists at `position` (another the anchor-split stage
/// iteration put one there — e.g. a shared cube corner anchor), we
/// reuse it without splitting any edge.
pub(super) fn insert_anchor_at_position(
    brep: &mut BRep,
    position: Point3,
    edge_splits: &mut Vec<(EdgeHandle, EdgeHandle, EdgeHandle)>,
) -> Result<VertexHandle, String> {
    for (vh, v) in brep.vertices.iter() {
        let d_sq = (v.position.x - position.x).powi(2)
            + (v.position.y - position.y).powi(2)
            + (v.position.z - position.z).powi(2);
        if d_sq < 1e-12 {
            return Ok(vh);
        }
    }
    if let Some((eh, t)) = find_edge_and_param_at(brep, position) {
        let split = brep
            .split_edge(eh, t)
            .map_err(|e| format!("split_edge for anchor failed: {e}"))?;
        edge_splits.push((eh, split.edge_a, split.edge_b));
        return Ok(split.vertex);
    }
    // Curve-aware fallback: the chord test above misses anchors on
    // CURVED adjacent edges twice over — perp-walk math offsets the
    // anchor along the tangent (off-curve by the chord sagitta), and
    // the containment test itself compares against the vertex chord
    // (measured: a boss-wall chain's anchors sat 5e-3 off the r=1
    // wall arc and the anchor-split stage failed "no edge contains anchor"). Find
    // the closest interior point on any curved edge within a
    // sagitta-scaled gate and split THERE — the realized vertex
    // lands exactly on the curve, which is more faithful than the
    // requested tangent-offset position.
    let mut best: Option<(EdgeHandle, Point3, f64)> = None;
    for (eh, edge) in brep.edges.iter() {
        let Some(curve) = edge.curve.as_ref() else {
            continue;
        };
        let (t0, t1) = (curve.t_start, curve.t_end);
        let mut bt = t0;
        let mut bd = f64::INFINITY;
        for k in 0..=32 {
            let t = t0 + (t1 - t0) * k as f64 / 32.0;
            let d = curve.geometry.point_at(t).distance_to(&position);
            if d < bd {
                bd = d;
                bt = t;
            }
        }
        let (mut lo, mut hi) = (
            (bt - (t1 - t0) / 32.0).max(t0),
            (bt + (t1 - t0) / 32.0).min(t1),
        );
        for _ in 0..30 {
            let m1 = lo + (hi - lo) / 3.0;
            let m2 = hi - (hi - lo) / 3.0;
            if curve.geometry.point_at(m1).distance_to(&position)
                < curve.geometry.point_at(m2).distance_to(&position)
            {
                hi = m2;
            } else {
                lo = m1;
            }
        }
        let t = 0.5 * (lo + hi);
        // Interior only — endpoint projections belong to the vertex
        // reuse path above.
        let frac = (t - t0) / (t1 - t0).max(1e-300);
        if !(0.001..=0.999).contains(&frac) {
            continue;
        }
        let proj = curve.geometry.point_at(t);
        let d = proj.distance_to(&position);
        if best.as_ref().is_none_or(|(_, _, bd2)| d < *bd2) {
            best = Some((eh, proj, d));
        }
    }
    let gate = brep.tolerance.linear().max(1.0e-9) * 1.0e4;
    if popkern_core::trace::trace_enabled("anchor-insert") {
        eprintln!(
            "[anchor-insert] position {position:?} best {:?} gate {gate:.3e}",
            best.as_ref().map(|(eh, p, d)| (*eh, (p.x, p.y, p.z), *d))
        );
    }
    if let Some((eh, proj, d)) = best
        && d < gate
    {
        let tol = brep.tolerance;
        let split = brep
            .split_edge_at_point(eh, proj, &tol)
            .map_err(|e| format!("curve-aware anchor split failed: {e}"))?;
        edge_splits.push((eh, split.edge_a, split.edge_b));
        return Ok(split.vertex);
    }
    // Loose vertex reuse: a NEIGHBOURING fillet's the anchor-split stage may have
    // already split the shared curved edge at (the on-curve image
    // of) this anchor — the perp-walk request sits a chord-sagitta
    // off the curve, so the tight reuse above misses it (measured on
    // the boss-wall chain: requested (0.9937, −0.0496) vs the
    // existing on-arc vertex 5.1e-3 away; the arc no longer has an
    // interior point to split). Reuse the nearest existing vertex
    // within the same sagitta-scaled gate.
    let mut best_v: Option<(VertexHandle, f64)> = None;
    for (vh, v) in brep.vertices.iter() {
        let d = v.position.distance_to(&position);
        if best_v.map(|(_, bd)| d < bd).unwrap_or(true) {
            best_v = Some((vh, d));
        }
    }
    if let Some((vh, d)) = best_v
        && d < gate
    {
        return Ok(vh);
    }
    Err(format!("no edge contains anchor position {position:?}"))
}

/// §28.H Session 3 (PK-141): land a genuinely INTERIOR anchor by
/// splicing an mev spur into the target face's loop from the nearest
/// boundary vertex. The spur edge carries the straight chord curve;
/// downstream phases treat the new vertex like any the anchor-split stage anchor
/// and the spur is consumed or dissolved by the corner machinery.
/// Transactional per spur: candidate clone, `validate()`-free (mid-
/// pipeline states may fail global checks) but structural mev
/// postconditions apply.
pub(in crate::construct) fn insert_anchor_with_spur(
    brep: &mut BRep,
    position: Point3,
    face: FaceHandle,
    tolerance: &Tolerance,
) -> Result<VertexHandle, String> {
    let gate_multiplier = 100.0;
    // The position must actually lie on the face's surface.
    let sof = brep
        .faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .ok_or_else(|| "spur target face has no surface".to_string())?;
    // Scoped to CURVED supports (PK-141): planar interior overrides
    // keep the documented Session-2 perp-walk fallback — the
    // direct-apex fixtures depend on it (measured: 9 corner tests
    // regressed when the spur landed planner apex targets).
    if matches!(sof.geometry, SurfaceGeometry::Plane { .. }) {
        return Err("interior anchor spur is scoped to curved supports".to_string());
    }
    use popkern_core::geometry::SurfaceEval;
    let on_surface =
        popkern_core::geometry::uv_projection::project_point_to_uv(&sof.geometry, &position)
            .map(|(u, v)| sof.geometry.point_at(u, v).distance_to(&position))
            .unwrap_or(f64::INFINITY);
    let gate = tolerance.linear().max(1e-9) * gate_multiplier;
    if on_surface > gate {
        return Err(format!(
            "anchor {position:?} is {on_surface:.3e} off face {face:?} — not an interior anchor"
        ));
    }
    // Nearest boundary half-edge origin across all loops.
    let fd = brep
        .faces
        .get(face)
        .ok_or_else(|| "spur target face vanished".to_string())?;
    let mut loops = vec![fd.outer_loop];
    loops.extend(fd.inner_loops.iter().copied());
    let mut best: Option<(popkern_core::topology::HalfEdgeHandle, f64)> = None;
    for lh in loops {
        let Some(l) = brep.loops.get(lh) else {
            continue;
        };
        let start = l.half_edge;
        let mut cur = start;
        for _ in 0..brep.half_edges.len() + 1 {
            let Some(he) = brep.half_edges.get(cur) else {
                break;
            };
            if let Some(v) = brep.vertices.get(he.origin) {
                let d = v.position.distance_to(&position);
                if best.map(|(_, bd)| d < bd).unwrap_or(true) {
                    best = Some((cur, d));
                }
            }
            cur = he.next;
            if cur == start {
                break;
            }
        }
    }
    let (he, _) = best.ok_or_else(|| "spur target face has no boundary vertices".to_string())?;
    let origin_pos = brep
        .half_edges
        .get(he)
        .and_then(|h| brep.vertices.get(h.origin))
        .map(|v| v.position)
        .ok_or_else(|| "spur origin vanished".to_string())?;
    let mev = brep
        .mev(he, position)
        .map_err(|e| format!("anchor spur mev failed: {e}"))?;
    if let Some(e) = brep.edges.get_mut(mev.edge) {
        e.curve = popkern_core::geometry::CurveOnEdge::line_segment(
            origin_pos,
            position,
            tolerance.linear(),
        );
    }
    if popkern_core::trace::trace_enabled("anchor") {
        eprintln!(
            "[anchor] interior anchor spur on {face:?}: {origin_pos:?} -> {position:?} (vertex {:?})",
            mev.vertex
        );
    }
    Ok(mev.vertex)
}

/// Splice an interior-rail spur at a KNOWN corner vertex
/// (PK-185/PK-186 reflex hole-corner ends). Unlike
/// [`insert_anchor_with_spur`] this is not scoped to curved
/// supports and does not search for the nearest boundary vertex:
/// the reflex classification already identified the exact corner
/// (the fillet edge's own vertex), and the spur must attach there —
/// it becomes the strip's end boundary through that corner. The
/// nearest-vertex search could attach to an unrelated feature
/// vertex, and the planar scope exists to protect direct-apex
/// planner targets, which never classify as interior-rail (their
/// edges run with `anchor_on_trim = false`).
pub(in crate::construct) fn insert_interior_rail_spur(
    brep: &mut BRep,
    face: FaceHandle,
    corner: VertexHandle,
    position: Point3,
    tolerance: &Tolerance,
) -> Result<VertexHandle, String> {
    let fd = brep
        .faces
        .get(face)
        .ok_or_else(|| format!("interior-rail spur target face {face:?} not found"))?;
    let mut loops = vec![fd.outer_loop];
    loops.extend(fd.inner_loops.iter().copied());
    let mut he_at_corner = None;
    'outer: for lh in loops {
        let Some(l) = brep.loops.get(lh) else {
            continue;
        };
        let start = l.half_edge;
        let mut cur = start;
        for _ in 0..brep.half_edges.len() + 1 {
            let Some(he) = brep.half_edges.get(cur) else {
                break;
            };
            if he.origin == corner {
                he_at_corner = Some(cur);
                break 'outer;
            }
            cur = he.next;
            if cur == start {
                break;
            }
        }
    }
    let he = he_at_corner.ok_or_else(|| {
        format!("interior-rail spur: corner {corner:?} not on any loop of face {face:?}")
    })?;
    let corner_pos = brep
        .vertices
        .get(corner)
        .map(|v| v.position)
        .ok_or_else(|| "interior-rail spur corner vanished".to_string())?;
    let mev = brep
        .mev(he, position)
        .map_err(|e| format!("interior-rail spur mev failed: {e}"))?;
    if let Some(e) = brep.edges.get_mut(mev.edge) {
        e.curve = popkern_core::geometry::CurveOnEdge::line_segment(
            corner_pos,
            position,
            tolerance.linear(),
        );
    }
    if popkern_core::trace::trace_enabled("anchor") {
        eprintln!(
            "[anchor] interior-rail spur on {face:?}: corner {corner:?} {corner_pos:?} -> {position:?} (vertex {:?})",
            mev.vertex
        );
    }
    Ok(mev.vertex)
}

/// Land one the anchor-split stage anchor. Non-reflex slots keep the historical
/// cascade (on-edge split → curved-support spur → perp-walk
/// fallback). Interior-rail slots first try an exact on-boundary
/// landing (when the rail end coincides with existing topology the
/// end is NOT reflex — the band reaches the boundary), then splice a
/// targeted spur at the recorded corner and register the reflex end
/// for the terminal-disc stage.
#[allow(clippy::too_many_arguments)]
pub(in crate::construct) fn land_anchor(
    brep: &mut BRep,
    requested: Point3,
    placement: super::perp_anchor_positions::AnchorPlacement,
    face: FaceHandle,
    corner: VertexHandle,
    at_start: bool,
    tolerance: &Tolerance,
    edge_splits: &mut Vec<(EdgeHandle, EdgeHandle, EdgeHandle)>,
    reflex_ends: &mut Vec<ReflexEndAnchor>,
) -> Result<VertexHandle, String> {
    if placement.junction {
        // The partner strip may have landed this section contact
        // already (vertex reuse); otherwise spur in from the junction.
        if let Ok(v) = insert_anchor_at_position(brep, requested, edge_splits) {
            return Ok(v);
        }
        return insert_interior_rail_spur(brep, face, corner, requested, tolerance);
    }
    if placement.interior_rail {
        if let Ok(v) = insert_anchor_at_position(brep, requested, edge_splits) {
            return Ok(v);
        }
        let v = insert_interior_rail_spur(brep, face, corner, requested, tolerance)?;
        reflex_ends.push(ReflexEndAnchor {
            face,
            at_start,
            corner,
            spur_tip: v,
        });
        return Ok(v);
    }
    insert_anchor_at_position(brep, requested, edge_splits)
        .or_else(|_| insert_anchor_with_spur(brep, requested, face, tolerance))
        .or_else(|_| insert_anchor_at_position(brep, placement.position, edge_splits))
}

/// A reflex hole-corner end landed with an interior-rail spur
/// (PK-185/PK-186). Consumed by the hole-rim terminal-disc stage
/// after the consumed-merge stage: the strip's end boundary currently runs
/// spur → corner → stub through the ORIGINAL sharp corner, and the
/// disc stage replaces the strip's share of that path with the
/// exact terminal cross-section arc, authoring the flat end face in
/// the transversal plane.
#[derive(Debug, Clone, Copy)]
pub struct ReflexEndAnchor {
    /// The support face that received the spur.
    pub face: FaceHandle,
    /// Whether this is the fillet edge's start vertex end.
    pub at_start: bool,
    /// The original sharp hole corner (the fillet edge's own vertex).
    pub corner: VertexHandle,
    /// The spur tip minted at the rail end (interior to `face`).
    pub spur_tip: VertexHandle,
}

/// Scan all edges in `brep` and find the one whose 3D line segment
/// contains `point` within a small tolerance, returning the edge
/// handle and parameter along the edge in `(0.001, 0.999)`.
///
/// `pub(crate)` for §6c-DirectApex.2.b — the apex builder calls this from
/// `resolve_apex_vertices`'s create branch when no existing BRep
/// vertex matches the planner's target position; the edge passing
/// through the target is then split via [`BRep::split_edge`] to
/// insert a new vertex on it.  Production-grade callers in this
/// crate (the anchor-split stage, multi_chord) already exercise this path on
/// cube fixtures.
pub(crate) fn find_edge_and_param_at(brep: &BRep, point: Point3) -> Option<(EdgeHandle, f64)> {
    for (eh, edata) in brep.edges.iter() {
        let he_h = edata.half_edge;
        let Some(he) = brep.half_edges.get(he_h) else {
            continue;
        };
        // The second vertex is the twin's origin on manifold
        // edges, or the next half-edge's origin for boundary
        // edges (whose twin is the null sentinel).
        let v1_pos = if let Some(twin) = brep.half_edges.get(he.twin) {
            brep.vertices[twin.origin].position
        } else if let Some(next_he) = brep.half_edges.get(he.next) {
            brep.vertices[next_he.origin].position
        } else {
            continue;
        };
        let v0 = brep.vertices[he.origin].position;
        let v1 = v1_pos;
        let dx = v1.x - v0.x;
        let dy = v1.y - v0.y;
        let dz = v1.z - v0.z;
        let len_sq = dx * dx + dy * dy + dz * dz;
        if len_sq < 1e-18 {
            continue;
        }
        let t = ((point.x - v0.x) * dx + (point.y - v0.y) * dy + (point.z - v0.z) * dz) / len_sq;
        if !(0.001..=0.999).contains(&t) {
            continue;
        }
        let px = v0.x + t * dx;
        let py = v0.y + t * dy;
        let pz = v0.z + t * dz;
        let dist_sq = (point.x - px).powi(2) + (point.y - py).powi(2) + (point.z - pz).powi(2);
        if dist_sq < 1e-10 {
            return Some((eh, t));
        }
    }
    None
}
