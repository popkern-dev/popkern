//! Boundary-completion ring rebuild for winding reconciliation
//! (PK-153): a face contorted across both winding classes keeps its
//! trusted live sides and closes its ring by authoring twin halves
//! against coincident dead sides lying on its own surface.

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{
    BRep, FaceHandle, HalfEdge, HalfEdgeHandle, LoopHandle, VertexHandle,
};

use super::{side_face, side_span, trace};

/// Rebuild a contorted face's outer ring by boundary completion:
/// order the trusted sides into chains, then close the ring by
/// authoring twin halves against coincident dead sides lying on the
/// face's carrier plane. The continuation walk is a bounded search
/// that accepts only a UNIQUE closure (PK-260): at a stall with
/// several coincident dead candidates the walk branches, and a ring
/// is adopted only when exactly one branch closes — correctness by
/// uniqueness, never by guessing. Returns seam vertices for
/// coincident-record merging, or None when the search stalls, fails
/// to close, or closes ambiguously.
pub(super) fn rebuild_face_ring(
    brep: &mut BRep,
    face: FaceHandle,
    trusted: &[HalfEdgeHandle],
) -> Option<Vec<VertexHandle>> {
    let tol = brep.tolerance.linear().max(1.0e-9) * 10.0;
    let loop_handle = brep.faces.get(face)?.outer_loop;
    // PK-260: surfaceless residual slivers (the demotion passes strip
    // their inherited carrier before rung 1 runs) still have a
    // well-defined planar carrier — their trusted sides' endpoints.
    // Without this the plane-gated span test rejects EVERY dead
    // continuation and sliver ring rebuilds always stall at zero
    // candidates (measured on the oblique base-corner class).
    let carrier = face_carrier_plane(brep, face, trusted, tol);
    let plan = plan_face_ring(brep, face, trusted, tol, carrier.as_ref())?;
    execute_face_ring(brep, face, loop_handle, &plan, tol)
}

/// One adopted continuation of the ring walk: an existing trusted
/// side, or a dead side whose face-side twin the execute phase
/// authors.
#[derive(Clone, Copy)]
enum RingStep {
    Trusted(HalfEdgeHandle),
    Dead(HalfEdgeHandle),
}

/// A unique closing walk found by [`plan_face_ring`], replayed
/// verbatim by [`execute_face_ring`].
struct RingPlan {
    first: HalfEdgeHandle,
    ring_start: VertexHandle,
    ring_start_pos: Point3,
    steps: Vec<RingStep>,
}

fn trusted_span(brep: &BRep, side: HalfEdgeHandle) -> Option<(VertexHandle, VertexHandle)> {
    let data = brep.half_edges.get(side)?;
    let twin = brep.half_edges.get(data.twin)?;
    Some((data.origin, twin.origin))
}

/// Pure search phase: walk trusted continuations greedily (they are
/// forced), branch at multi-candidate dead-side stalls, and accept
/// only a unique closing walk. No topology is mutated here, so
/// branches need no rollback.
fn plan_face_ring(
    brep: &BRep,
    face: FaceHandle,
    trusted: &[HalfEdgeHandle],
    tol: f64,
    carrier: Option<&(Point3, popkern_core::geometry::Vec3)>,
) -> Option<RingPlan> {
    let first = trusted.first().copied()?;
    let (ring_start, first_end) = trusted_span(brep, first)?;
    let ring_start_pos = brep.vertices.get(ring_start)?.position;
    let position = |vertex: VertexHandle| brep.vertices.get(vertex).map(|v| v.position);

    struct Node {
        steps: Vec<RingStep>,
        pool: Vec<HalfEdgeHandle>,
        used_dead: Vec<HalfEdgeHandle>,
        open_end: VertexHandle,
    }
    let depth_cap = brep.half_edges.len() + trusted.len() + 8;
    let mut stack = vec![Node {
        steps: Vec::new(),
        pool: trusted[1..].to_vec(),
        used_dead: Vec::new(),
        open_end: first_end,
    }];
    let mut solution: Option<Vec<RingStep>> = None;
    let mut visits = 0usize;
    // Deepest stall (position, candidate count, depth) for the
    // decline diagnostics.
    let mut deepest_stall: Option<(Point3, usize, usize)> = None;
    while let Some(mut node) = stack.pop() {
        visits += 1;
        if visits > 512 {
            trace(&format!(
                "winding reconciliation: ring rebuild for {face:?} exceeded its search budget"
            ));
            return None;
        }
        loop {
            if node.steps.len() + 1 > depth_cap {
                break;
            }
            let Some(open_pos) = position(node.open_end) else {
                break;
            };
            if open_pos.distance_to(&ring_start_pos) <= tol && node.steps.len() + 1 > 2 {
                if solution.is_some() {
                    trace(&format!(
                        "winding reconciliation: ring rebuild for {face:?} closed ambiguously"
                    ));
                    return None;
                }
                solution = Some(node.steps.clone());
                break;
            }
            // Trusted continuation first (forced, as in the original
            // greedy walk).
            if let Some(index) = node.pool.iter().position(|side| {
                trusted_span(brep, *side).is_some_and(|(start, _)| {
                    position(start).is_some_and(|p| p.distance_to(&open_pos) <= tol)
                })
            }) {
                let side = node.pool.remove(index);
                let Some((_, end)) = trusted_span(brep, side) else {
                    break;
                };
                node.steps.push(RingStep::Trusted(side));
                node.open_end = end;
                continue;
            }
            // Dead-side continuations on the carrier.
            let candidates: Vec<HalfEdgeHandle> = brep
                .half_edges
                .iter()
                .filter(|(handle, data)| {
                    !brep.half_edges.contains_key(data.twin) && !node.used_dead.contains(handle)
                })
                .filter(|(handle, _)| side_face(brep, *handle) != Some(face))
                .filter(|(handle, _)| {
                    let Some((start, end)) = side_span(brep, *handle) else {
                        return false;
                    };
                    end.distance_to(&open_pos) <= tol
                        && span_on_carrier_plane(carrier, &start, &end, tol)
                })
                .map(|(handle, _)| handle)
                .collect();
            match candidates.as_slice() {
                [] => {
                    let depth = node.steps.len();
                    if deepest_stall.is_none_or(|(_, _, d)| depth > d) {
                        deepest_stall = Some((open_pos, 0, depth));
                    }
                    break;
                }
                [only] => {
                    let Some(origin) = brep.half_edges.get(*only).map(|data| data.origin) else {
                        break;
                    };
                    node.used_dead.push(*only);
                    node.steps.push(RingStep::Dead(*only));
                    node.open_end = origin;
                    continue;
                }
                multiple => {
                    for side in multiple {
                        let Some(origin) = brep.half_edges.get(*side).map(|data| data.origin)
                        else {
                            continue;
                        };
                        let mut steps = node.steps.clone();
                        steps.push(RingStep::Dead(*side));
                        let mut used_dead = node.used_dead.clone();
                        used_dead.push(*side);
                        stack.push(Node {
                            steps,
                            pool: node.pool.clone(),
                            used_dead,
                            open_end: origin,
                        });
                    }
                    break;
                }
            }
        }
    }
    match solution {
        Some(steps) => Some(RingPlan {
            first,
            ring_start,
            ring_start_pos,
            steps,
        }),
        None => {
            if let Some((pos, count, _)) = deepest_stall {
                trace(&format!(
                    "winding reconciliation: ring rebuild for {face:?} found {count} dead \
                     continuation(s) at ({:.3},{:.3},{:.3})",
                    pos.x, pos.y, pos.z
                ));
                // PK-260 stall census: every side incident to the
                // stall position, with liveness/ownership — the
                // measured basis for widening the continuation
                // sourcing.
                if crate::direct_apex_trace::enabled() {
                    for (handle, data) in brep.half_edges.iter() {
                        let Some((start, end)) = side_span(brep, handle) else {
                            continue;
                        };
                        if start.distance_to(&pos) > tol && end.distance_to(&pos) > tol {
                            continue;
                        }
                        let twin_live = brep.half_edges.contains_key(data.twin);
                        let owner = side_face(brep, handle);
                        trace(&format!(
                            "  stall side {handle:?} face={owner:?} twin_live={twin_live} \
                             span=({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3})",
                            start.x, start.y, start.z, end.x, end.y, end.z
                        ));
                    }
                }
            } else {
                trace(&format!(
                    "winding reconciliation: ring rebuild for {face:?} did not close"
                ));
            }
            None
        }
    }
}

/// Mutating replay of a unique ring plan: author the face-side twins
/// for adopted dead sides, collect seam vertices, and wire the chain.
fn execute_face_ring(
    brep: &mut BRep,
    face: FaceHandle,
    loop_handle: LoopHandle,
    plan: &RingPlan,
    tol: f64,
) -> Option<Vec<VertexHandle>> {
    let mut chain = vec![plan.first];
    let (_, mut open_end) = trusted_span(brep, plan.first)?;
    let mut seam_vertices = Vec::new();
    for step in &plan.steps {
        match step {
            RingStep::Trusted(side) => {
                let (start, end) = trusted_span(brep, *side)?;
                if start != open_end {
                    seam_vertices.push(start);
                    seam_vertices.push(open_end);
                }
                chain.push(*side);
                open_end = end;
            }
            RingStep::Dead(dead) => {
                let dead_data = brep.half_edges.get(*dead)?;
                let dead_edge = dead_data.edge;
                let dead_origin = dead_data.origin;
                let dead_end = brep
                    .half_edges
                    .get(dead_data.next)
                    .map(|next| next.origin)?;
                let new_side = brep.half_edges.insert(HalfEdge {
                    origin: dead_end,
                    twin: *dead,
                    next: Default::default(),
                    prev: Default::default(),
                    edge: dead_edge,
                    loop_: loop_handle,
                    pcurve: None,
                });
                brep.half_edges[*dead].twin = new_side;
                brep.boundary_edges.remove(&dead_edge);
                if dead_end != open_end {
                    seam_vertices.push(dead_end);
                    seam_vertices.push(open_end);
                }
                chain.push(new_side);
                open_end = dead_origin;
            }
        }
    }
    // The plan guarantees closure; re-verify cheaply so a stale
    // handle between plan and execute rolls back via the caller.
    let open_pos = brep.vertices.get(open_end)?.position;
    if open_pos.distance_to(&plan.ring_start_pos) > tol || chain.len() < 3 {
        trace(&format!(
            "winding reconciliation: ring rebuild for {face:?} did not close ({} side(s))",
            chain.len()
        ));
        return None;
    }
    if open_end != plan.ring_start {
        seam_vertices.push(open_end);
        seam_vertices.push(plan.ring_start);
    }
    let n = chain.len();
    for k in 0..n {
        let side = chain[k];
        let next = chain[(k + 1) % n];
        brep.half_edges[side].next = next;
        brep.half_edges[next].prev = side;
        brep.half_edges[side].loop_ = loop_handle;
        brep.half_edges[side].pcurve = None;
    }
    brep.loops[loop_handle].half_edge = chain[0];
    Some(seam_vertices)
}

/// The face's planar carrier: its own Plane surface when present,
/// otherwise (surfaceless residual sliver) the plane fitted through
/// its trusted sides' endpoints — accepted only when every endpoint
/// actually lies on that plane, so the walk stays conservative.
/// Non-planar surfaces yield None (dead-side completion declines).
fn face_carrier_plane(
    brep: &BRep,
    face: FaceHandle,
    trusted: &[HalfEdgeHandle],
    tol: f64,
) -> Option<(Point3, popkern_core::geometry::Vec3)> {
    if let Some(sof) = brep.faces.get(face).and_then(|d| d.surface.as_ref()) {
        let SurfaceGeometry::Plane { origin, normal, .. } = &sof.geometry else {
            return None;
        };
        return normal.normalized().map(|n| (*origin, n));
    }
    let mut points = Vec::new();
    for side in trusted {
        let Some(data) = brep.half_edges.get(*side) else {
            continue;
        };
        if let Some(v) = brep.vertices.get(data.origin) {
            points.push(v.position);
        }
        if let Some(twin) = brep.half_edges.get(data.twin)
            && let Some(v) = brep.vertices.get(twin.origin)
        {
            points.push(v.position);
        }
    }
    if points.len() < 3 {
        return None;
    }
    // Newell's method over the point cloud.
    let mut n = popkern_core::geometry::Vec3::new(0.0, 0.0, 0.0);
    for i in 0..points.len() {
        let a = points[i];
        let b = points[(i + 1) % points.len()];
        n.x += (a.y - b.y) * (a.z + b.z);
        n.y += (a.z - b.z) * (a.x + b.x);
        n.z += (a.x - b.x) * (a.y + b.y);
    }
    let n = n.normalized()?;
    let origin = points[0];
    let coplanar = points.iter().all(|p| ((*p - origin).dot(&n)).abs() <= tol);
    if !coplanar {
        return None;
    }
    Some((origin, n))
}

/// Both span endpoints lie on the carrier plane.
fn span_on_carrier_plane(
    carrier: Option<&(Point3, popkern_core::geometry::Vec3)>,
    start: &Point3,
    end: &Point3,
    tol: f64,
) -> bool {
    let Some((origin, n)) = carrier else {
        return false;
    };
    let dist = |p: &Point3| ((*p - *origin).dot(n)).abs();
    dist(start) <= tol && dist(end) <= tol
}

/// Chord-ring completion (PK-214): after the ghost complex retires
/// and the mirror class reverses, the rewritten support quads still
/// carry PARTNERLESS dead diagonal chords (the shortcut sides the
/// contorted interface authored) while their true far paths lie as
/// dead sides on the neighbor faces. For each planar face whose dead
/// sides include a partnerless chord, sever the chords and close the
/// ring by boundary completion against the coincident dead sides on
/// its own carrier. Returns (rebuilt_rings, seam_vertices).
pub(in super::super) fn rebuild_partnerless_chord_rings(
    brep: &mut BRep,
) -> (usize, Vec<VertexHandle>) {
    let tol = brep.tolerance.linear().max(1.0e-9) * 10.0;
    let has_coincident_counterpart = |brep: &BRep, side: HalfEdgeHandle| -> bool {
        let Some((a, b)) = side_span(brep, side) else {
            return false;
        };
        brep.half_edges.iter().any(|(other, data)| {
            if other == side || brep.half_edges.contains_key(data.twin) {
                return false;
            }
            let Some((c, d)) = side_span(brep, other) else {
                return false;
            };
            (a.distance_to(&c) <= tol && b.distance_to(&d) <= tol)
                || (a.distance_to(&d) <= tol && b.distance_to(&c) <= tol)
        })
    };
    let mut rebuilt = 0usize;
    let mut seam_vertices = Vec::new();
    let faces: Vec<FaceHandle> = brep.faces.keys().collect();
    for face in faces {
        let sides = super::face_sides(brep, face);
        if sides.is_empty() {
            continue;
        }
        let chords: Vec<HalfEdgeHandle> = sides
            .iter()
            .copied()
            .filter(|side| {
                brep.half_edges
                    .get(*side)
                    .is_some_and(|data| !brep.half_edges.contains_key(data.twin))
                    && !has_coincident_counterpart(brep, *side)
            })
            .collect();
        if chords.is_empty() {
            continue;
        }
        let trusted: Vec<HalfEdgeHandle> = sides
            .iter()
            .copied()
            .filter(|side| {
                brep.half_edges
                    .get(*side)
                    .is_some_and(|data| brep.half_edges.contains_key(data.twin))
            })
            .collect();
        if trusted.is_empty() {
            trace(&format!(
                "chord-ring completion declined for {face:?}: no trusted live sides"
            ));
            continue;
        }
        let snapshot = brep.clone();
        for chord in &chords {
            super::remove_side(brep, *chord);
        }
        match rebuild_face_ring(brep, face, &trusted) {
            Some(mut seam) => {
                rebuilt += 1;
                seam_vertices.append(&mut seam);
            }
            None => {
                trace(&format!(
                    "chord-ring completion: ring rebuild declined for {face:?}, rolled back"
                ));
                *brep = snapshot;
            }
        }
    }
    (rebuilt, seam_vertices)
}
