//! Carrier absorption of coplanar residual scaffolding.
//!
//! After the V-free partition, planar production faces (base, facets)
//! are separated from their true seam partners by retired scaffold
//! fragments that lie ON the same plane carrier: the band between a
//! facet's old strip boundary and the new generator, and the 2-gon
//! slivers holding displaced base spans. Per the OCC reference the
//! production faces extend across those regions. This pass absorbs
//! each such fragment into its production neighbor with the
//! transactional `merge_faces_across_edge` (dV=0 dE=-1 dF=-1 dL=-1),
//! iterated to fixpoint.
//!
//! Guards: the fragment's sides must all fit the production plane
//! (endpoints and curve midpoint within the linear tolerance), and
//! every fragment vertex must already appear on some production loop —
//! fragments still wired to scaffold-interior topology (e.g. the
//! retired apex vertex spurs) are left for the apex-cell teardown.

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::TopologyRole;

use super::*;

/// Absorb coplanar residual fragments into adjacent planar production
/// faces until no more qualify; returns how many rewrites were applied.
pub(super) fn absorb_coplanar_residual_scaffold(brep: &mut BRep) -> usize {
    let mut rewrites = 0usize;
    let mut declined: std::collections::BTreeSet<HalfEdgeHandle> =
        std::collections::BTreeSet::new();
    // Retired scaffolding may already fail whole-model validation for
    // reasons this pass neither causes nor can repair; the acceptance
    // criterion for each rewrite is "no NEW validation failure".
    let preexisting = brep.validate().err().map(|err| err.to_string());
    if let Some(reason) = preexisting.as_deref()
        && crate::direct_apex_trace::enabled()
    {
        eprintln!(
            "[direct-apex-trace §14.9/carrier-absorb] pre-existing validation failure \
             (rewrites gated on not changing it): {reason}"
        );
    }
    let acceptable = |outcome: Result<(), String>| -> Result<(), String> {
        match (outcome, preexisting.as_deref()) {
            (Ok(()), _) => Ok(()),
            (Err(err), Some(prior)) if err == prior => Ok(()),
            (Err(err), _) => Err(err),
        }
    };
    loop {
        // Bigon welds first: they don't need coplanarity (the fragment's
        // two sides are geometrically the same span) and they unblock
        // merge candidates by removing degenerate slivers.
        let candidate = next_bigon_weld_side(brep, &declined);
        if let Some(production_side) = candidate {
            let snapshot = brep.clone();
            let far_twin_loop_closed_before = bigon_weld_operands(brep, production_side)
                .and_then(|operands| {
                    let loop_ = brep.half_edges.get(operands.far_twin)?.loop_;
                    Some(collect_loop_half_edges(brep, loop_).is_ok())
                })
                .unwrap_or(false);
            let outcome = weld_across_bigon(brep, production_side).and_then(|()| {
                local_weld_postconditions(brep, production_side, far_twin_loop_closed_before)
                    .map_err(|err| format!("post-weld check: {err}"))
            });
            match outcome {
                Ok(()) => {
                    rewrites += 1;
                    continue;
                }
                Err(reason) => {
                    if crate::direct_apex_trace::enabled() {
                        eprintln!(
                            "[direct-apex-trace §14.9/carrier-absorb] bigon weld at \
                             {production_side:?} rolled back: {reason}"
                        );
                        if let Some(operands) = bigon_weld_operands(&snapshot, production_side) {
                            eprintln!(
                                "[direct-apex-trace §14.9/carrier-absorb]   operands: near \
                                 {:?} far {:?} far_twin {:?} (far_twin loop {:?})",
                                operands.near_side,
                                operands.far_side,
                                operands.far_twin,
                                snapshot
                                    .half_edges
                                    .get(operands.far_twin)
                                    .map(|data| data.loop_),
                            );
                        }
                    }
                    *brep = snapshot;
                    declined.insert(production_side);
                    continue;
                }
            }
        }
        let Some(merge_side) = next_absorbable_side(brep, &declined) else {
            break;
        };
        let snapshot = brep.clone();
        let outcome = brep
            .merge_faces_across_edge(merge_side)
            .map_err(|err| err.to_string())
            .and_then(|_| {
                acceptable(brep.validate().map_err(|err| err.to_string()))
                    .map_err(|err| format!("post-merge validation: {err}"))
            });
        match outcome {
            Ok(()) => {
                rewrites += 1;
                if crate::direct_apex_trace::enabled() {
                    eprintln!(
                        "[direct-apex-trace §14.9/carrier-absorb] merged residual fragment \
                         across {merge_side:?}"
                    );
                }
            }
            Err(err) => {
                *brep = snapshot;
                declined.insert(merge_side);
                if crate::direct_apex_trace::enabled() {
                    eprintln!(
                        "[direct-apex-trace §14.9/carrier-absorb] merge across {merge_side:?} \
                         rolled back: {err}"
                    );
                }
            }
        }
    }
    rewrites
}

/// A production half-edge whose twin lives in a residual two-sided
/// fragment (a bigon sliver holding a displaced duplicate of the span).
fn next_bigon_weld_side(
    brep: &BRep,
    declined: &std::collections::BTreeSet<HalfEdgeHandle>,
) -> Option<HalfEdgeHandle> {
    for (_, face_data) in &brep.faces {
        if face_data.role == TopologyRole::ConstructionResidual {
            continue;
        }
        let Ok(cycle) = collect_loop_half_edges(brep, face_data.outer_loop) else {
            continue;
        };
        for side in cycle {
            if !declined.contains(&side) && bigon_weld_operands(brep, side).is_some() {
                return Some(side);
            }
        }
    }
    None
}

struct BigonWeldOperands {
    fragment_face: FaceHandle,
    fragment_loop: popkern_core::topology::LoopHandle,
    near_side: HalfEdgeHandle,
    far_side: HalfEdgeHandle,
    far_twin: HalfEdgeHandle,
}

fn bigon_weld_operands(brep: &BRep, production_side: HalfEdgeHandle) -> Option<BigonWeldOperands> {
    let data = brep.half_edges.get(production_side)?;
    let near_side = data.twin;
    let near = brep.half_edges.get(near_side)?;
    let fragment_loop = near.loop_;
    let fragment_face = brep.loops.get(fragment_loop)?.face;
    let fragment = brep.faces.get(fragment_face)?;
    if fragment.role != TopologyRole::ConstructionResidual || !fragment.inner_loops.is_empty() {
        return None;
    }
    let cycle = collect_loop_half_edges(brep, fragment_loop).ok()?;
    if cycle.len() != 2 || !cycle.contains(&near_side) {
        return None;
    }
    let far_side = cycle.into_iter().find(|he| *he != near_side)?;
    let far = brep.half_edges.get(far_side)?;
    let far_twin = far.twin;
    // Scaffolding carries non-mutual twin pointers (measured: the far
    // side pointing back into the bigon). Only a mutual partner outside
    // the fragment can become the production twin.
    if far_twin == near_side || far_twin == far_side || far_twin == production_side {
        return None;
    }
    let far_twin_data = brep.half_edges.get(far_twin)?;
    if far_twin_data.twin != far_side {
        return None;
    }
    // Never adopt an unwired boundary record as a production twin.
    brep.loops.get(far_twin_data.loop_)?;
    // A twin must traverse opposite to the production side: it starts
    // where the production side ends — the near side's origin.
    (far_twin_data.origin == near.origin).then_some(())?;
    // Displaced-membership guard: scaffold chains from OTHER loops may
    // route through the bigon's half-edges (loop_ pointers and
    // next-cycles disagree on retired records). Removing the pair would
    // break those chains, so decline until the deeper handover retires
    // them.
    let removed = [near_side, far_side];
    for (other, other_data) in &brep.half_edges {
        if removed.contains(&other) {
            continue;
        }
        if removed.contains(&other_data.next) || removed.contains(&other_data.prev) {
            return None;
        }
    }
    for (other_loop, loop_data) in &brep.loops {
        if other_loop != fragment_loop && removed.contains(&loop_data.half_edge) {
            return None;
        }
    }
    Some(BigonWeldOperands {
        fragment_face,
        fragment_loop,
        near_side,
        far_side,
        far_twin,
    })
}

/// Retire a residual bigon sliver by re-pairing the production side
/// with the far side's twin: `dHE=-2 dE=-1 dF=-1 dL=-1`, Euler-sound.
fn weld_across_bigon(brep: &mut BRep, production_side: HalfEdgeHandle) -> Result<(), String> {
    let operands = bigon_weld_operands(brep, production_side)
        .ok_or_else(|| "bigon operands vanished".to_owned())?;
    let production_edge = brep
        .half_edges
        .get(production_side)
        .map(|data| data.edge)
        .ok_or_else(|| "stale production side".to_owned())?;
    let far_edge = brep
        .half_edges
        .get(operands.far_side)
        .map(|data| data.edge)
        .ok_or_else(|| "stale far side".to_owned())?;
    // Re-pair.
    brep.half_edges[production_side].twin = operands.far_twin;
    brep.half_edges[operands.far_twin].twin = production_side;
    brep.half_edges[operands.far_twin].edge = production_edge;
    if brep
        .edges
        .get(production_edge)
        .is_some_and(|edge| edge.half_edge == operands.near_side)
    {
        brep.edges[production_edge].half_edge = production_side;
    }
    // Retire the fragment and its now-orphaned records.
    let removed = [operands.near_side, operands.far_side];
    for stale in removed {
        brep.half_edges.remove(stale);
    }
    repair_vertex_outgoing_pointers(brep, production_side, operands.far_twin);
    if far_edge != production_edge {
        brep.edges.remove(far_edge);
        brep.boundary_edges.remove(&far_edge);
    }
    // Displaced scaffold members may still claim the fragment loop
    // (loop_ pointers and anchor cycles disagree on retired records):
    // re-anchor the loop to a surviving claimer instead of orphaning
    // them; delete the loop and face only when nothing claims them.
    let survivor = brep
        .half_edges
        .iter()
        .find_map(|(heh, data)| (data.loop_ == operands.fragment_loop).then_some(heh));
    if let Some(anchor) = survivor {
        brep.loops[operands.fragment_loop].half_edge = anchor;
    } else {
        brep.loops.remove(operands.fragment_loop);
        brep.faces.remove(operands.fragment_face);
    }
    if crate::direct_apex_trace::enabled() {
        eprintln!(
            "[direct-apex-trace §14.9/carrier-absorb] bigon {:?} welded: {production_side:?} \
             re-paired with {:?}",
            operands.fragment_face, operands.far_twin
        );
    }
    Ok(())
}

/// Local structural postconditions for a bigon weld: the re-paired
/// sides are mutual twins on one edge and both owning loops still close
/// through them.
fn local_weld_postconditions(
    brep: &BRep,
    production_side: HalfEdgeHandle,
    far_twin_loop_closed_before: bool,
) -> Result<(), String> {
    let data = brep
        .half_edges
        .get(production_side)
        .ok_or_else(|| "production side vanished".to_owned())?;
    let twin = brep
        .half_edges
        .get(data.twin)
        .ok_or_else(|| "new twin vanished".to_owned())?;
    if twin.twin != production_side {
        return Err("twin pairing is not mutual".to_owned());
    }
    if twin.edge != data.edge {
        return Err("twin pair spans two edge records".to_owned());
    }
    // The production loop must stay intact; the residual partner loop
    // only needs to not get WORSE (scaffold chains may be pre-broken).
    for (side, must_close) in [
        (production_side, true),
        (data.twin, far_twin_loop_closed_before),
    ] {
        if !must_close {
            continue;
        }
        let loop_handle = brep
            .half_edges
            .get(side)
            .map(|d| d.loop_)
            .ok_or_else(|| "side vanished".to_owned())?;
        let cycle = collect_loop_half_edges(brep, loop_handle)
            .map_err(|_| format!("loop {loop_handle:?} no longer closes"))?;
        if !cycle.contains(&side) {
            return Err(format!(
                "loop {loop_handle:?} dropped {side:?} from its cycle"
            ));
        }
    }
    Ok(())
}

/// Any vertex whose `outgoing` pointer died with the bigon gets a live
/// replacement whose origin matches.
fn repair_vertex_outgoing_pointers(
    brep: &mut BRep,
    production_side: HalfEdgeHandle,
    far_twin: HalfEdgeHandle,
) {
    let candidates = [production_side, far_twin];
    let vertex_handles: Vec<VertexHandle> = brep.vertices.keys().collect();
    for vertex in vertex_handles {
        let outgoing = brep.vertices[vertex].outgoing;
        if brep.half_edges.contains_key(outgoing) {
            continue;
        }
        if let Some(replacement) = candidates
            .iter()
            .copied()
            .find(|he| brep.half_edges.get(*he).is_some_and(|d| d.origin == vertex))
            .or_else(|| {
                brep.half_edges
                    .iter()
                    .find_map(|(heh, data)| (data.origin == vertex).then_some(heh))
            })
        {
            brep.vertices[vertex].outgoing = replacement;
        }
    }
}

/// A production-plane half-edge whose twin's face is an absorbable
/// coplanar residual fragment.
fn next_absorbable_side(
    brep: &BRep,
    declined: &std::collections::BTreeSet<HalfEdgeHandle>,
) -> Option<HalfEdgeHandle> {
    let production_vertices = production_loop_vertices(brep);
    for (face, face_data) in &brep.faces {
        if face_data.role == TopologyRole::ConstructionResidual {
            continue;
        }
        let Some(surface) = face_data.surface.as_ref() else {
            continue;
        };
        if !matches!(surface.geometry, SurfaceGeometry::Plane { .. }) {
            continue;
        }
        let Ok(cycle) = collect_loop_half_edges(brep, face_data.outer_loop) else {
            continue;
        };
        for side in cycle {
            if declined.contains(&side) {
                continue;
            }
            let Some(twin) = brep
                .half_edges
                .get(side)
                .and_then(|data| brep.half_edges.get(data.twin))
            else {
                continue;
            };
            let Some(neighbor_loop) = brep.loops.get(twin.loop_) else {
                continue;
            };
            let neighbor = neighbor_loop.face;
            if neighbor == face {
                continue;
            }
            let Some(neighbor_data) = brep.faces.get(neighbor) else {
                continue;
            };
            if neighbor_data.role != TopologyRole::ConstructionResidual
                || !neighbor_data.inner_loops.is_empty()
            {
                continue;
            }
            match fragment_absorbability(brep, neighbor, &surface.geometry, &production_vertices) {
                Ok(()) => return Some(side),
                Err(reason) => {
                    if crate::direct_apex_trace::enabled() {
                        eprintln!(
                            "[direct-apex-trace §14.9/carrier-absorb] fragment {neighbor:?} \
                             (next to {face:?}) not absorbable: {reason}"
                        );
                    }
                }
            }
        }
    }
    None
}

fn production_loop_vertices(brep: &BRep) -> std::collections::BTreeSet<VertexHandle> {
    let mut vertices = std::collections::BTreeSet::new();
    for (_, face_data) in &brep.faces {
        if face_data.role == TopologyRole::ConstructionResidual {
            continue;
        }
        let Ok(cycle) = collect_loop_half_edges(brep, face_data.outer_loop) else {
            continue;
        };
        for side in cycle {
            if let Some(data) = brep.half_edges.get(side) {
                vertices.insert(data.origin);
            }
        }
    }
    vertices
}

/// Every fragment side fits the plane carrier and every fragment vertex
/// already participates in production topology.
fn fragment_absorbability(
    brep: &BRep,
    fragment: FaceHandle,
    carrier: &SurfaceGeometry,
    production_vertices: &std::collections::BTreeSet<VertexHandle>,
) -> Result<(), String> {
    let fragment_data = brep
        .faces
        .get(fragment)
        .ok_or_else(|| "stale fragment".to_owned())?;
    let cycle = collect_loop_half_edges(brep, fragment_data.outer_loop)
        .map_err(|_| "fragment loop does not close".to_owned())?;
    let fit_tolerance = brep.tolerance.linear().max(1.0e-9) * 10.0;
    for side in cycle {
        let data = brep
            .half_edges
            .get(side)
            .ok_or_else(|| format!("stale side {side:?}"))?;
        if !production_vertices.contains(&data.origin) {
            return Err(format!(
                "vertex {:?} is scaffold-interior (not on any production loop)",
                data.origin
            ));
        }
        let origin = brep
            .vertices
            .get(data.origin)
            .ok_or_else(|| format!("stale vertex on {side:?}"))?;
        let vertex_gap = point_plane_distance(carrier, origin.position);
        if vertex_gap > fit_tolerance {
            return Err(format!(
                "vertex {:?} off the carrier by {vertex_gap:.3e}",
                data.origin
            ));
        }
        if let Some(curve) = brep
            .edges
            .get(data.edge)
            .and_then(|edge| edge.curve.as_ref())
        {
            let midpoint = curve.geometry.point_at(0.5 * (curve.t_start + curve.t_end));
            let midpoint_gap = point_plane_distance(carrier, midpoint);
            if midpoint_gap > fit_tolerance {
                return Err(format!(
                    "edge {:?} midpoint off the carrier by {midpoint_gap:.3e}",
                    data.edge
                ));
            }
        }
    }
    Ok(())
}

fn point_plane_distance(carrier: &SurfaceGeometry, point: popkern_core::geometry::Point3) -> f64 {
    let Some((u, v)) = carrier.closest_point(&point) else {
        return f64::INFINITY;
    };
    carrier.point_at(u, v).distance_to(&point)
}
