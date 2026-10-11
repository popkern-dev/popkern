//! Redundant shared-boundary CHAINS between co-surfacial faces
//! (PK-379).
//!
//! [`super::redundant_edge`] dissolves an edge only when it is the sole
//! edge its two faces share: dissolving one edge of a multi-edge
//! shared boundary leaves the survivors traversed twice by the merged
//! loop (a zero-width slit). When the shared edges form ONE PATH whose
//! interior vertices carry nothing but the path, the whole boundary is
//! dissolvable soundly: `kef` the first edge, then every remaining
//! path edge is a spur from the freed end and retires with `kev`. The
//! merged loop is exactly `A ∪ B` minus the path.
//!
//! Measured: the two-axis staggered cube union kept its coplanar top
//! and bottom skins in three pieces each (A-only, overlap, B-only)
//! because the overlap rectangle shares two edges meeting at an
//! interior corner with each neighbour; the x-only stagger, whose
//! pieces share single edges, merged to six faces.

use popkern_core::geometry::intersection::coincident_analytic::are_surfaces_coincident_analytic;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};
use std::collections::BTreeMap;

/// The two faces an edge separates, when they are distinct,
/// co-surfacial and equally oriented (outward normals, not flags).
fn mergeable_face_pair(
    brep: &BRep,
    eh: EdgeHandle,
    tolerance: &Tolerance,
) -> Option<(FaceHandle, FaceHandle)> {
    let edge = brep.edges.get(eh)?;
    let he = brep.half_edges.get(edge.half_edge)?;
    if he.twin == edge.half_edge {
        return None;
    }
    let twin = brep.half_edges.get(he.twin)?;
    let face_a = brep.loops.get(he.loop_)?.face;
    let face_b = brep.loops.get(twin.loop_)?.face;
    if face_a == face_b {
        return None;
    }
    let surf_a = brep.faces.get(face_a)?.surface.as_ref()?;
    let surf_b = brep.faces.get(face_b)?.surface.as_ref()?;
    if !are_surfaces_coincident_analytic(&surf_a.geometry, &surf_b.geometry, tolerance)
        || !super::redundant_edge::outward_normals_agree(brep, eh, surf_a, surf_b)?
    {
        return None;
    }
    Some((face_a, face_b))
}

fn endpoints(brep: &BRep, eh: EdgeHandle) -> Option<(VertexHandle, VertexHandle)> {
    let edge = brep.edges.get(eh)?;
    let he = brep.half_edges.get(edge.half_edge)?;
    let twin = brep.half_edges.get(he.twin)?;
    Some((he.origin, twin.origin))
}

fn valence(brep: &BRep, v: VertexHandle) -> usize {
    brep.half_edges
        .iter()
        .filter(|(_, he)| he.origin == v)
        .count()
}

/// The ordered path of edges `eh`'s two co-surfacial faces share, when
/// that shared boundary is one path (two ends, interior vertices used
/// by nothing else). A single shared edge is a path of length one.
pub(crate) fn shared_coplanar_chain(
    brep: &BRep,
    eh: EdgeHandle,
    tolerance: &Tolerance,
) -> Option<Vec<EdgeHandle>> {
    let (face_a, face_b) = mergeable_face_pair(brep, eh, tolerance)?;
    let shared: Vec<EdgeHandle> = brep
        .edges
        .iter()
        .filter_map(|(oh, other)| {
            let ohe = brep.half_edges.get(other.half_edge)?;
            let otwin = brep.half_edges.get(ohe.twin)?;
            let fa = brep.loops.get(ohe.loop_)?.face;
            let fb = brep.loops.get(otwin.loop_)?.face;
            ((fa == face_a && fb == face_b) || (fa == face_b && fb == face_a)).then_some(oh)
        })
        .collect();
    if shared.len() == 1 {
        return Some(shared);
    }
    // Degree census over the shared edges.
    let mut degree: BTreeMap<VertexHandle, Vec<EdgeHandle>> = BTreeMap::new();
    for &e in &shared {
        let (a, b) = endpoints(brep, e)?;
        degree.entry(a).or_default().push(e);
        degree.entry(b).or_default().push(e);
    }
    let ends: Vec<VertexHandle> = degree
        .iter()
        .filter(|(_, es)| es.len() == 1)
        .map(|(v, _)| *v)
        .collect();
    if ends.len() != 2 || degree.values().any(|es| es.len() > 2) {
        return None;
    }
    // Interior vertices must carry only the two path edges, or the
    // merged loop would still pass through them.
    for (v, es) in &degree {
        if es.len() == 2 && valence(brep, *v) != 2 {
            return None;
        }
    }
    // Order the path from the smaller end handle (deterministic).
    let mut path = Vec::with_capacity(shared.len());
    let mut at = ends[0].min(ends[1]);
    let mut prev: Option<EdgeHandle> = None;
    while path.len() < shared.len() {
        let next = degree
            .get(&at)?
            .iter()
            .copied()
            .find(|e| Some(*e) != prev)?;
        let (a, b) = endpoints(brep, next)?;
        at = if a == at { b } else { a };
        prev = Some(next);
        path.push(next);
    }
    Some(path)
}

/// Dissolve one shared chain: `kef` the first edge, then `kev` each
/// remaining edge from its freed (valence-1) end; the merged face's
/// pcurves are re-derived (the killed face's half-edges carried
/// pcurves in ITS chart — coincident planes need not share an origin
/// or u-axis, and a stale chart strands the next boolean's records;
/// measured on the flush-arm gusset union → lip-slot subtract).
/// Returns the surviving face. Fails when the chain does not collapse
/// as predicted; the caller stages.
pub(crate) fn dissolve_chain(brep: &mut BRep, chain: &[EdgeHandle]) -> Result<FaceHandle, String> {
    let Some((&first, rest)) = chain.split_first() else {
        return Err("empty chain".into());
    };
    let he = brep.edges.get(first).ok_or("stale chain edge")?.half_edge;
    let kef = brep.kef(he).map_err(|e| format!("kef: {e}"))?;
    let survivor = brep
        .loops
        .get(kef.surviving_loop)
        .ok_or("stale surviving loop")?
        .face;
    for &e in rest {
        let (a, b) = endpoints(brep, e).ok_or("stale chain edge")?;
        let edge = brep.edges.get(e).ok_or("stale chain edge")?;
        let he = edge.half_edge;
        let twin = brep.half_edges.get(he).ok_or("stale half-edge")?.twin;
        // kev kills the vertex at twin.origin — pick the half-edge
        // arriving at the valence-1 end.
        let arriving = if valence(brep, b) == 1 {
            he
        } else if valence(brep, a) == 1 {
            twin
        } else {
            return Err(format!("chain edge {e:?} has no free end after kef"));
        };
        brep.kev(arriving).map_err(|err| format!("kev: {err}"))?;
    }
    reset_face_pcurves(brep, survivor)?;
    Ok(survivor)
}

/// Clear every pcurve on `face`'s loops and re-derive them in the
/// face's own chart.
fn reset_face_pcurves(brep: &mut BRep, face: FaceHandle) -> Result<(), String> {
    let f = brep.faces.get(face).ok_or("stale merged face")?;
    let mut loops = vec![f.outer_loop];
    loops.extend(f.inner_loops.iter().copied());
    let mut hes = Vec::new();
    for lh in loops {
        let start = brep.loops.get(lh).ok_or("stale loop")?.half_edge;
        let mut cursor = start;
        for _ in 0..=brep.half_edges.len() {
            hes.push(cursor);
            cursor = brep.half_edges.get(cursor).ok_or("stale half-edge")?.next;
            if cursor == start {
                break;
            }
        }
    }
    for h in hes {
        if let Some(he) = brep.half_edges.get_mut(h) {
            he.pcurve = None;
        }
    }
    brep.recompute_face_pcurves(face)
        .map_err(|e| format!("pcurve recompute on merged face: {e}"))
}
