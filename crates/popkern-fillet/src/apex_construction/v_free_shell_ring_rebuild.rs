//! Boundary-completion ring rebuild for winding reconciliation
//! (PK-153): a face contorted across both winding classes keeps its
//! trusted live sides and closes its ring by authoring twin halves
//! against coincident dead sides lying on its own surface.

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, FaceHandle, HalfEdge, HalfEdgeHandle, VertexHandle};

use super::{side_face, side_span, trace};

/// Rebuild a contorted face's outer ring by boundary completion:
/// order the trusted sides into chains, then close the ring by
/// authoring twin halves against coincident dead sides lying on the
/// face's own surface. Returns seam vertices for coincident-record
/// merging, or None when the walk is ambiguous or fails to close.
pub(super) fn rebuild_face_ring(
    brep: &mut BRep,
    face: FaceHandle,
    trusted: &[HalfEdgeHandle],
) -> Option<Vec<VertexHandle>> {
    let tol = brep.tolerance.linear().max(1.0e-9) * 10.0;
    let loop_handle = brep.faces.get(face)?.outer_loop;
    let span_of = |brep: &BRep, side: HalfEdgeHandle| -> Option<(VertexHandle, VertexHandle)> {
        let data = brep.half_edges.get(side)?;
        let twin = brep.half_edges.get(data.twin)?;
        Some((data.origin, twin.origin))
    };
    let mut chain: Vec<HalfEdgeHandle> = Vec::new();
    let mut pool: Vec<HalfEdgeHandle> = trusted.to_vec();
    let first = pool.first().copied()?;
    pool.remove(0);
    chain.push(first);
    let (ring_start, mut open_end) = span_of(brep, first)?;
    let ring_start_pos = brep.vertices.get(ring_start)?.position;
    let mut seam_vertices = Vec::new();
    for _ in 0..brep.half_edges.len() + trusted.len() + 8 {
        let open_pos = brep.vertices.get(open_end)?.position;
        if open_pos.distance_to(&ring_start_pos) <= tol && chain.len() > 2 {
            break;
        }
        // Trusted continuation first.
        if let Some(index) = pool.iter().position(|side| {
            span_of(brep, *side).is_some_and(|(start, _)| {
                brep.vertices
                    .get(start)
                    .is_some_and(|v| v.position.distance_to(&open_pos) <= tol)
            })
        }) {
            let side = pool.remove(index);
            let (start, end) = span_of(brep, side)?;
            if start != open_end {
                seam_vertices.push(start);
                seam_vertices.push(open_end);
            }
            chain.push(side);
            open_end = end;
            continue;
        }
        // Dead-side completion: a dead side ending at the open
        // position whose span lies on this face's surface gets its
        // face-side twin authored here.
        let candidates: Vec<HalfEdgeHandle> = brep
            .half_edges
            .iter()
            .filter(|(_, data)| !brep.half_edges.contains_key(data.twin))
            .filter(|(handle, _)| side_face(brep, *handle) != Some(face))
            .filter(|(handle, _)| {
                let Some((start, end)) = side_span(brep, *handle) else {
                    return false;
                };
                end.distance_to(&open_pos) <= tol
                    && span_on_face_surface(brep, face, &start, &end, tol)
            })
            .map(|(handle, _)| handle)
            .collect();
        let chosen = match candidates.as_slice() {
            [only] => Some(*only),
            [] => None,
            multiple => {
                // PK-260: coincident duplicate dead sides (the
                // far-side trim pair — the strip's stolen-trim slot
                // and the support replacement's slot) both end at the
                // open position. Resolve by ring-closure lookahead
                // instead of declining: prefer the candidate whose far
                // end closes the ring or chains onward (trusted pool
                // first, then another dead side); tie-break to the
                // side whose edge carries an exact curve (the authored
                // trim) over a curve-less span. A residual tie still
                // declines — the walk must not guess.
                let far_pos = |side: HalfEdgeHandle| -> Option<Point3> {
                    let origin = brep.half_edges.get(side)?.origin;
                    brep.vertices.get(origin).map(|v| v.position)
                };
                let score = |side: HalfEdgeHandle| -> Option<usize> {
                    let pos = far_pos(side)?;
                    let closes = pos.distance_to(&ring_start_pos) <= tol;
                    let trusted_next = pool.iter().any(|other| {
                        span_of(brep, *other).is_some_and(|(start, _)| {
                            brep.vertices
                                .get(start)
                                .is_some_and(|v| v.position.distance_to(&pos) <= tol)
                        })
                    });
                    let dead_next = brep.half_edges.iter().any(|(handle, data)| {
                        handle != side
                            && !brep.half_edges.contains_key(data.twin)
                            && side_span(brep, handle)
                                .is_some_and(|(_, end)| end.distance_to(&pos) <= tol)
                    });
                    let curve = brep
                        .half_edges
                        .get(side)
                        .and_then(|data| brep.edges.get(data.edge))
                        .is_some_and(|edge| edge.curve.is_some());
                    let rank = if closes {
                        3
                    } else if trusted_next {
                        2
                    } else if dead_next {
                        1
                    } else {
                        0
                    };
                    Some(rank * 2 + usize::from(curve))
                };
                let mut best: Option<(usize, HalfEdgeHandle)> = None;
                let mut tied = false;
                for side in multiple {
                    let Some(s) = score(*side) else { continue };
                    match &best {
                        Some((b, _)) if s == *b => tied = true,
                        Some((b, _)) if s > *b => {
                            best = Some((s, *side));
                            tied = false;
                        }
                        None => best = Some((s, *side)),
                        _ => {}
                    }
                }
                match (best, tied) {
                    (Some((_, side)), false) => Some(side),
                    _ => None,
                }
            }
        };
        let Some(dead) = chosen else {
            trace(&format!(
                "winding reconciliation: ring rebuild for {face:?} found {} dead \
                 continuation(s) at ({:.3},{:.3},{:.3})",
                candidates.len(),
                open_pos.x,
                open_pos.y,
                open_pos.z
            ));
            return None;
        };
        let dead_data = brep.half_edges.get(dead)?;
        let dead_edge = dead_data.edge;
        let dead_origin = dead_data.origin;
        let dead_end = brep
            .half_edges
            .get(dead_data.next)
            .map(|next| next.origin)?;
        let new_side = brep.half_edges.insert(HalfEdge {
            origin: dead_end,
            twin: dead,
            next: Default::default(),
            prev: Default::default(),
            edge: dead_edge,
            loop_: loop_handle,
            pcurve: None,
        });
        brep.half_edges[dead].twin = new_side;
        brep.boundary_edges.remove(&dead_edge);
        if dead_end != open_end {
            seam_vertices.push(dead_end);
            seam_vertices.push(open_end);
        }
        chain.push(new_side);
        open_end = dead_origin;
    }
    let open_pos = brep.vertices.get(open_end)?.position;
    if open_pos.distance_to(&ring_start_pos) > tol || chain.len() < 3 {
        trace(&format!(
            "winding reconciliation: ring rebuild for {face:?} did not close ({} side(s))",
            chain.len()
        ));
        return None;
    }
    if open_end != ring_start {
        seam_vertices.push(open_end);
        seam_vertices.push(ring_start);
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

/// Both span endpoints lie on the face's surface (plane distance for
/// planes; other geometries decline so the walk stays conservative).
fn span_on_face_surface(
    brep: &BRep,
    face: FaceHandle,
    start: &Point3,
    end: &Point3,
    tol: f64,
) -> bool {
    let Some(sof) = brep.faces.get(face).and_then(|d| d.surface.as_ref()) else {
        return false;
    };
    let SurfaceGeometry::Plane { origin, normal, .. } = &sof.geometry else {
        return false;
    };
    let Some(n) = normal.normalized() else {
        return false;
    };
    let dist = |p: &Point3| ((*p - *origin).dot(&n)).abs();
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
