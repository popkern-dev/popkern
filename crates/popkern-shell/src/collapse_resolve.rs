//! Face-collapse detection and Euler dissolution (PK-326 row 477).
//!
//! A face DIES at an offset distance when no sample of its offset
//! surface realizes the offset distance against the original body's
//! boundary — every sample sits closer to some OTHER boundary face
//! than to its own generator, so the Minkowski result contains no
//! trace of the face (a thin rib under erosion past half its width).
//! Detection is measured (trimmed distance queries, no flag algebra);
//! resolution dissolves each dead component through Euler operators
//! (`kef` merges, `kev` spur retirement, terminal `kvr`) on a
//! candidate skin, committing only when the component reduces
//! completely and the survivor validates. Anything unresolved is left
//! for the chord-inversion collapse gate to decline — folded geometry
//! never serves.

use std::collections::{BTreeMap, BTreeSet, HashMap};

use crate::body_distance::BodyDistanceIndex;
use crate::error::{ShellError, ShellWarning};
use crate::skin_map::SkinMap;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

/// Grid resolution for per-face death sampling.
const GRID: usize = 5;

/// Detect dead faces and dissolve their components. Returns the
/// `FaceCollapsed` warnings for resolved components; unresolvable
/// components leave the skin untouched (the downstream chord gate
/// declines them).
///
/// # Errors
///
/// Only on internal inconsistencies ([`ShellError::Topology`]);
/// geometric infeasibility is not an error here.
pub(crate) fn resolve_face_collapse(
    outer: &BRep,
    skin: &mut BRep,
    skin_map: &SkinMap,
    offset_surfaces: &HashMap<FaceHandle, SurfaceGeometry>,
    tolerance: &Tolerance,
) -> Result<Vec<ShellWarning>, ShellError> {
    let index = BodyDistanceIndex::build(outer);
    // Per-face death tests are pure reads of the index (GRID² samples ×
    // min-over-all-faces distance queries — the stage's cost driver);
    // dispatch in parallel with an ordered collect (PK-329) so the dead
    // set keeps the index's sorted face order.
    let verdicts = crate::parallel_dispatch::dispatch_ordered(&index.entries, false, |entry| {
        let Some(offset) = offset_surfaces.get(&entry.face) else {
            return false;
        };
        face_is_dead(entry, offset, &index, tolerance)
    });
    let dead_outer: Vec<FaceHandle> = index
        .entries
        .iter()
        .zip(verdicts)
        .filter_map(|(entry, dead)| dead.then_some(entry.face))
        .collect();
    if dead_outer.is_empty() {
        return Ok(Vec::new());
    }

    // Skin handles of the dead set.
    let mut dead_skin: BTreeSet<FaceHandle> = BTreeSet::new();
    let mut skin_to_outer: BTreeMap<FaceHandle, FaceHandle> = BTreeMap::new();
    for ofh in &dead_outer {
        if let Some(&sfh) = skin_map.faces.get(ofh)
            && skin.faces.contains_key(sfh)
        {
            dead_skin.insert(sfh);
            skin_to_outer.insert(sfh, *ofh);
        }
    }

    let mut warnings = Vec::new();
    for component in components(skin, &dead_skin) {
        let mut candidate = skin.clone();
        match dissolve_component(&mut candidate, &component) {
            Ok(()) => {
                let mut faces: Vec<FaceHandle> = component
                    .iter()
                    .filter_map(|sfh| skin_to_outer.get(sfh).copied())
                    .collect();
                faces.sort();
                warnings.push(ShellWarning::FaceCollapsed {
                    faces,
                    mechanism: "no offset sample realizes the offset distance against \
                                the body boundary (feature thinner than the offset); \
                                dissolved via kef/kev/kvr into the surviving skin"
                        .to_string(),
                });
                *skin = candidate;
            }
            Err(_mechanism) => {
                // Unresolved: leave the skin alone; the chord-inversion
                // gate declines with OffsetTooLarge as before.
            }
        }
    }
    Ok(warnings)
}

/// A face is dead when at least one grid sample is valid and none
/// realizes its own offset distance against the body boundary.
fn face_is_dead(
    entry: &crate::body_distance::FaceTrimEntry,
    offset: &SurfaceGeometry,
    index: &BodyDistanceIndex,
    tolerance: &Tolerance,
) -> bool {
    let Some((u0, v0, u1, v1)) = entry.uv_bbox() else {
        return false;
    };
    let mut valid = 0usize;
    for i in 0..GRID {
        for j in 0..GRID {
            let u = u0 + (u1 - u0) * (i as f64 + 0.5) / GRID as f64;
            let v = v0 + (v1 - v0) * (j as f64 + 0.5) / GRID as f64;
            if !entry.contains(u, v) {
                continue;
            }
            let q0 = entry.surface.point_at(u, v);
            let Some(sample) = crate::curved_recompute::project_to_surface(offset, &q0) else {
                continue;
            };
            let d_f = sample.distance_to(&q0);
            if d_f <= tolerance.linear() * 10.0 {
                return false; // zero-offset face: trivially alive
            }
            let Some(dist) = index.distance(&sample) else {
                continue;
            };
            valid += 1;
            let band = (tolerance.linear() * 100.0).max(d_f * 1.0e-3);
            if dist >= d_f - band {
                return false; // this sample survives: the face is alive
            }
        }
    }
    valid > 0
}

/// Connected components of the dead set under shared edges.
fn components(skin: &BRep, dead: &BTreeSet<FaceHandle>) -> Vec<Vec<FaceHandle>> {
    let mut adjacency: BTreeMap<FaceHandle, BTreeSet<FaceHandle>> = BTreeMap::new();
    for (_, edge) in &skin.edges {
        let Some(he) = skin.half_edges.get(edge.half_edge) else {
            continue;
        };
        let Some(twin) = skin.half_edges.get(he.twin) else {
            continue;
        };
        let (Some(l1), Some(l2)) = (skin.loops.get(he.loop_), skin.loops.get(twin.loop_)) else {
            continue;
        };
        if dead.contains(&l1.face) && dead.contains(&l2.face) && l1.face != l2.face {
            adjacency.entry(l1.face).or_default().insert(l2.face);
            adjacency.entry(l2.face).or_default().insert(l1.face);
        }
    }
    let mut seen: BTreeSet<FaceHandle> = BTreeSet::new();
    let mut out = Vec::new();
    for &f in dead {
        if seen.contains(&f) {
            continue;
        }
        let mut stack = vec![f];
        let mut comp = Vec::new();
        while let Some(g) = stack.pop() {
            if !seen.insert(g) {
                continue;
            }
            comp.push(g);
            if let Some(nbs) = adjacency.get(&g) {
                stack.extend(nbs.iter().copied());
            }
        }
        comp.sort();
        out.push(comp);
    }
    out
}

/// Dissolve one dead component on the candidate skin: merge interior
/// dead-dead edges (`kef`), retire seam spurs (`kev`), then absorb the
/// final plug face into its single alive neighbor and retire the ring.
fn dissolve_component(skin: &mut BRep, component: &[FaceHandle]) -> Result<(), String> {
    let mut dead: BTreeSet<FaceHandle> = component.iter().copied().collect();
    for f in &dead {
        let Some(face) = skin.faces.get(*f) else {
            return Err("dead face missing from skin".into());
        };
        if !face.inner_loops.is_empty() {
            return Err("dead face carries inner loops (pierced feature)".into());
        }
    }

    // Reduce the component to a single plug face.
    let bound = skin.edges.len() * 2 + 8;
    for _ in 0..bound {
        if let Some(he) = find_dead_dead_edge(skin, &dead) {
            let killed = skin
                .half_edges
                .get(skin.half_edges[he].twin)
                .and_then(|t| skin.loops.get(t.loop_))
                .map(|l| l.face)
                .ok_or("dead-dead edge lost its twin face")?;
            skin.kef(he).map_err(|e| format!("collapse kef: {e}"))?;
            dead.remove(&killed);
        } else if let Some(he) = find_seam_spur(skin, &dead) {
            skin.kev(he).map_err(|e| format!("collapse kev: {e}"))?;
        } else {
            break;
        }
    }
    if dead.len() != 1 {
        return Err(format!(
            "component did not reduce to one plug face ({} remain)",
            dead.len()
        ));
    }
    let plug = *dead.iter().next().ok_or("empty component")?;

    // Every remaining plug edge must border ONE alive face, and the
    // shared boundary must be an inner loop of that face (a blister
    // across several faces needs re-splitting — PK-327 territory).
    let (ring_he, alive, alive_loop) = plug_ring(skin, plug)?;
    let alive_face = skin
        .faces
        .get(alive)
        .ok_or("plug neighbor missing from skin")?;
    if !alive_face.inner_loops.contains(&alive_loop) {
        return Err("plug ring is not an inner loop of its alive neighbor".into());
    }

    // Absorb the plug (alive side survives), then retire the slit.
    let anchor = skin.half_edges[ring_he].origin;
    skin.kef(ring_he).map_err(|e| format!("plug kef: {e}"))?;
    for _ in 0..bound {
        if let Some(he) = find_loop_spur(skin, alive_loop) {
            skin.kev(he).map_err(|e| format!("slit kev: {e}"))?;
        } else {
            break;
        }
    }
    // Terminal state: empty ring + orphan anchor vertex.
    let loop_data = skin
        .loops
        .get(alive_loop)
        .ok_or("plug ring loop vanished before kvr")?;
    if skin.half_edges.contains_key(loop_data.half_edge) {
        return Err("plug slit did not retire completely".into());
    }
    skin.kvr(alive, alive_loop, anchor)
        .map_err(|e| format!("terminal kvr: {e}"))?;
    Ok(())
}

/// A half-edge separating two DISTINCT dead faces (sorted scan).
fn find_dead_dead_edge(
    skin: &BRep,
    dead: &BTreeSet<FaceHandle>,
) -> Option<popkern_core::topology::HalfEdgeHandle> {
    let mut edges: Vec<EdgeHandle> = skin.edges.keys().collect();
    edges.sort();
    for eh in edges {
        let he = skin.edges[eh].half_edge;
        let he_data = skin.half_edges.get(he)?;
        let twin = skin.half_edges.get(he_data.twin)?;
        let f1 = skin.loops.get(he_data.loop_)?.face;
        let f2 = skin.loops.get(twin.loop_)?.face;
        if f1 != f2 && dead.contains(&f1) && dead.contains(&f2) {
            return Some(he);
        }
    }
    None
}

/// A retire-able seam spur inside a dead face: a half-edge immediately
/// followed by its twin whose tip vertex has valence 1.
fn find_seam_spur(
    skin: &BRep,
    dead: &BTreeSet<FaceHandle>,
) -> Option<popkern_core::topology::HalfEdgeHandle> {
    let mut edges: Vec<EdgeHandle> = skin.edges.keys().collect();
    edges.sort();
    for eh in edges {
        let he = skin.edges[eh].half_edge;
        let he_data = skin.half_edges.get(he)?;
        let twin_h = he_data.twin;
        let twin = skin.half_edges.get(twin_h)?;
        let f1 = skin.loops.get(he_data.loop_)?.face;
        let f2 = skin.loops.get(twin.loop_)?.face;
        if f1 != f2 || !dead.contains(&f1) {
            continue;
        }
        for (spur, spur_twin) in [(he, twin_h), (twin_h, he)] {
            if skin.half_edges[spur].next == spur_twin {
                let tip = skin.half_edges[spur_twin].origin;
                let valence = skin.half_edges.values().filter(|h| h.origin == tip).count();
                if valence == 1 {
                    return Some(spur);
                }
            }
        }
    }
    None
}

/// A spur half-edge inside the given loop (slit retirement).
fn find_loop_spur(
    skin: &BRep,
    loop_: popkern_core::topology::LoopHandle,
) -> Option<popkern_core::topology::HalfEdgeHandle> {
    let loop_data = skin.loops.get(loop_)?;
    let start = loop_data.half_edge;
    if !skin.half_edges.contains_key(start) {
        return None;
    }
    let mut cursor = start;
    for _ in 0..skin.half_edges.len() {
        let he = skin.half_edges.get(cursor)?;
        if he.next == he.twin {
            let tip = skin.half_edges[he.twin].origin;
            let valence = skin.half_edges.values().filter(|h| h.origin == tip).count();
            if valence == 1 {
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

/// The plug's single ring: verify every plug edge borders the same
/// alive face and return (a ring half-edge on the ALIVE side, the
/// alive face, the alive-side loop).
fn plug_ring(
    skin: &BRep,
    plug: FaceHandle,
) -> Result<
    (
        popkern_core::topology::HalfEdgeHandle,
        FaceHandle,
        popkern_core::topology::LoopHandle,
    ),
    String,
> {
    let face = skin.faces.get(plug).ok_or("plug face missing")?;
    let start = skin
        .loops
        .get(face.outer_loop)
        .map(|l| l.half_edge)
        .ok_or("plug loop missing")?;
    let mut alive: Option<(popkern_core::topology::HalfEdgeHandle, FaceHandle, _)> = None;
    let mut cursor = start;
    for _ in 0..skin.half_edges.len() {
        let he = skin.half_edges.get(cursor).ok_or("plug loop broken")?;
        let twin = skin.half_edges.get(he.twin).ok_or("plug twin missing")?;
        let twin_loop = twin.loop_;
        let nb = skin
            .loops
            .get(twin_loop)
            .map(|l| l.face)
            .ok_or("plug neighbor loop missing")?;
        if nb != plug {
            match alive {
                None => alive = Some((he.twin, nb, twin_loop)),
                Some((_, f, _)) if f == nb => {}
                Some(_) => {
                    return Err("plug ring spans multiple alive faces (blister across \
                                edges; re-splitting is the trim stage's job)"
                        .into());
                }
            }
        }
        cursor = he.next;
        if cursor == start {
            break;
        }
    }
    alive.ok_or_else(|| "plug has no alive neighbor (erosion consumed the whole skin)".into())
}
