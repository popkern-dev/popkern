//! Hole-filling faces dissolved into their host (plan 4.8, PK-984).
//!
//! A flush union can leave a face whose whole boundary is one inner
//! loop of a co-surfacial neighbour: a block's top with a pocket
//! opening, and the plug's top filling it exactly. The two lie on one
//! oriented surface and the ring between them bounds nothing — the
//! redundant-edge pass leaves it (a single shared edge merges by `kef`,
//! a shared open path by chain dissolution; a closed ring would need
//! the host's inner loop removed, not a merge).
//!
//! The ring dissolves by Euler operators only, the sequence feature
//! dissolution uses: `kef` on one ring edge absorbs the plug into the
//! host (the plug's loop splices into the host's hole loop as a slit
//! path traversed both ways), `kev` retires that path from its leaf
//! ends, and `kvr` removes the empty ring with its last vertex. Each
//! dissolution is staged on a copy and kept only when the body still
//! validates; a ring vertex that also carries another edge (a wall
//! edge meeting the rim) stops the `kev` retirement short and the copy
//! is dropped.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle};

use crate::serve_region::ServeRegion;

/// The plug's ring on its host: (a ring half-edge on the host side, the
/// host face, the host's inner loop), when every edge of `plug`'s only
/// loop borders the same inner loop of one other face on the same
/// oriented surface.
fn plug_on_host(
    brep: &BRep,
    plug: FaceHandle,
    tolerance: &Tolerance,
) -> Option<(HalfEdgeHandle, FaceHandle, LoopHandle)> {
    use popkern_core::geometry::intersection::coincident_analytic::are_surfaces_coincident_analytic;
    let face = brep.faces.get(plug)?;
    if !face.inner_loops.is_empty() {
        return None;
    }
    let start = brep.loops.get(face.outer_loop)?.half_edge;
    let mut host: Option<(HalfEdgeHandle, LoopHandle)> = None;
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        let he = brep.half_edges.get(cursor)?;
        let twin = he.mate().and_then(|m| brep.half_edges.get(m))?;
        if he.mate() == Some(cursor) {
            return None;
        }
        match host {
            None => host = Some((he.mate()?, twin.loop_)),
            Some((_, l)) if l == twin.loop_ => {}
            Some(_) => return None,
        }
        cursor = he.next;
        if cursor == start {
            break;
        }
    }
    let (ring_he, ring_loop) = host?;
    let host_face = brep.loops.get(ring_loop)?.face;
    let host_data = brep.faces.get(host_face)?;
    if host_face == plug || !host_data.inner_loops.contains(&ring_loop) {
        return None;
    }
    let (sa, sb) = (host_data.surface.as_ref()?, face.surface.as_ref()?);
    if !are_surfaces_coincident_analytic(&sa.geometry, &sb.geometry, tolerance)
        || !super::merge_rechart::same_chart(&sa.geometry, &sb.geometry, tolerance)
    {
        return None;
    }
    let edge = brep.half_edges.get(ring_he)?.edge;
    if !super::redundant_edge::outward_normals_agree(brep, edge, sa, sb)? {
        return None;
    }
    Some((ring_he, host_face, ring_loop))
}

/// A half-edge of `loop_` whose successor is its own twin, ending at a
/// vertex no other half-edge leaves: the tip of the slit path.
fn slit_tip(brep: &BRep, loop_: LoopHandle) -> Option<HalfEdgeHandle> {
    let start = brep.loops.get(loop_)?.half_edge;
    if !brep.half_edges.contains_key(start) {
        return None;
    }
    let mut cursor = start;
    for _ in 0..brep.half_edges.len() {
        let he = brep.half_edges.get(cursor)?;
        if he.mate() == Some(he.next) {
            let tip = he.mate().and_then(|m| brep.half_edges.get(m))?.origin;
            if brep.half_edges.values().filter(|h| h.origin == tip).count() == 1 {
                return Some(cursor);
            }
        }
        cursor = he.next;
        if cursor == start {
            break;
        }
    }
    None
}

/// Absorb `plug` into its host along the ring and retire the ring.
fn dissolve(
    brep: &mut BRep,
    ring_he: HalfEdgeHandle,
    host: FaceHandle,
    ring: LoopHandle,
) -> Option<()> {
    let anchor = brep.half_edges.get(ring_he)?.origin;
    brep.kef(ring_he).ok()?;
    for _ in 0..=brep.edges.len() {
        match slit_tip(brep, ring) {
            Some(he) => {
                brep.kev(he).ok()?;
            }
            None => break,
        }
    }
    brep.kvr(host, ring, anchor).ok()
}

/// Dissolve every hole-filling face `region` admits into its host.
/// Returns how many it dissolved.
pub(crate) fn dissolve_filled_holes(
    result: &mut BRep,
    tolerance: &Tolerance,
    region: &ServeRegion,
) -> usize {
    let mut merged = 0usize;
    loop {
        let candidates: Vec<FaceHandle> = result.faces.keys().collect();
        let mut progressed = false;
        for plug in candidates {
            let Some((ring_he, host, ring)) = plug_on_host(result, plug, tolerance) else {
                continue;
            };
            let edge = result.half_edges[ring_he].edge;
            if !region.admits_edge(result, edge) {
                continue;
            }
            let mut staged = result.clone();
            if dissolve(&mut staged, ring_he, host, ring).is_some() && staged.validate().is_ok() {
                *result = staged;
                merged += 1;
                progressed = true;
                break;
            }
        }
        if !progressed {
            break;
        }
    }
    merged
}

#[cfg(test)]
#[path = "redundant_ring_tests.rs"]
mod tests;
