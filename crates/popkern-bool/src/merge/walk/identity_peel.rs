//! Peeled holes: a planar face's far, untouched holes skip its walk
//! (PK-802).
//!
//! A 64-slot disc's cap faces carry 63 slot holes; each later cut
//! chart-walked all of them again (~3 ms per cap, every cut) to
//! rediscover loops no record touches. The walk instead sees the
//! face's remaining half-edges, and each peeled hole goes back to the
//! one fragment that can enclose it ([`reattach_peeled`]); a face where
//! that is not decided re-walks whole, as without peeling.
//!
//! A peeled fragment carries no chart view (its peeled holes were not
//! chart-walked); chart views are read only by the torus sampler, and
//! only planar faces peel.

use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::topology::{BRep, FaceHandle, LoopHandle};

use super::super::graph::{EdgeSource, FaceRef, IntersectionGraph, Operand};
use super::identity::{OperandHalfEdgeMaps, closes_in_graph, mapped_loop};
use super::{FaceWalk, Fragment, walk_face};

/// Walk one face, peeling its far untouched holes when `peel` holds
/// them; `compare` (`POPKERN_BOOL_SCOPE=compare`) also walks the whole
/// face, serves that walk and asserts the peeled result matches it.
pub(super) fn walk_face_peeling(
    ctx: &mut FaceWalk<'_, '_>,
    face_ref: &FaceRef,
    he_indices: &[usize],
    peel: Option<PeeledHoles>,
    compare: bool,
) {
    let Some(peel) = peel else {
        walk_face(ctx, face_ref, he_indices);
        return;
    };
    let (start, fallbacks) = (ctx.fragments.len(), ctx.walk_fallbacks.len());
    walk_face(ctx, face_ref, &peel.rest);
    let reattached = reattach_peeled(ctx.graph, &mut ctx.fragments[start..], &peel, ctx.tolerance);
    if reattached && !compare {
        return;
    }
    // Not decided (or compare mode): walk the whole face.
    let peeled: Vec<Fragment> = ctx.fragments.drain(start..).collect();
    ctx.walk_fallbacks.truncate(fallbacks);
    ctx.chart_views.remove(face_ref);
    walk_face(ctx, face_ref, he_indices);
    if reattached {
        compare_peeled(&peeled, &ctx.fragments[start..]);
    }
}

/// Inner loops of a planar face peeled off before its walk (see
/// [`peel_far_holes`]).
pub(super) struct PeeledHoles {
    /// The peeled loops, each as the walk would seed it.
    pub loops: Vec<Vec<usize>>,
    /// The face's remaining half-edges, in their graph order.
    pub rest: Vec<usize>,
}

/// Peel the untouched holes of a PLANAR face that lie away from the
/// opposite operand, so the walk sees only the loops a cut can change.
///
/// A 64-slot disc's top face carries 63 slot holes; each later cut
/// re-walked all of them (~3 ms per cap face, every cut). A hole is
/// peeled when its loop is untouched (as for [`identity_fragment`]),
/// shares no node with the face's other half-edges (the walk cannot
/// route into or out of it), and its conservative curve box is
/// disjoint from `opposite_box` padded by 100× the linear tolerance:
/// on a plane, a loop enclosed by another lies in that loop's box, and
/// every record loop lies in the opposite operand's box, so no record
/// loop can enclose a peeled hole. Whether the hole belongs to the one
/// fragment the reduced walk yields is decided by
/// [`reattach_peeled`]; any other outcome re-walks the whole face.
#[allow(clippy::too_many_arguments)]
pub(super) fn peel_far_holes(
    graph: &IntersectionGraph,
    face_ref: &FaceRef,
    he_indices: &[usize],
    brep_a: &BRep,
    brep_b: &BRep,
    maps: &OperandHalfEdgeMaps<'_>,
    overlap_pairs: &[(FaceHandle, FaceHandle)],
    opposite_box: Option<Aabb>,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Option<PeeledHoles> {
    use popkern_core::geometry::surface::SurfaceGeometry;
    let opposite_box = opposite_box?;
    let (brep, map) = match face_ref.operand {
        Operand::A => (brep_a, maps.a),
        Operand::B => (brep_b, maps.b),
    };
    let face = brep.faces.get(face_ref.face)?;
    if face.inner_loops.is_empty()
        || !matches!(
            face.surface.as_ref().map(|s| &s.geometry),
            Some(SurfaceGeometry::Plane { .. })
        )
        || overlap_pairs
            .iter()
            .any(|&(fa, fb)| match face_ref.operand {
                Operand::A => fa == face_ref.face,
                Operand::B => fb == face_ref.face,
            })
        || graph.seam_substitutes.iter().any(|(f, _)| f == face_ref)
    {
        return None;
    }
    let on_face: std::collections::HashSet<usize> = he_indices.iter().copied().collect();
    let pad = 100.0 * tolerance.linear();
    let mut candidates: Vec<(LoopHandle, Vec<usize>)> = Vec::new();
    for &lh in &face.inner_loops {
        let Some(lp) = mapped_loop(brep, lh, map) else {
            continue;
        };
        let untouched = lp.iter().all(|&h| {
            let he = &graph.half_edges[h];
            on_face.contains(&h)
                && he.face.as_ref() == Some(face_ref)
                && matches!(
                    he.source,
                    EdgeSource::OriginalA(_) | EdgeSource::OriginalB(_)
                )
        }) && closes_in_graph(graph, &lp);
        if untouched && loop_box(brep, lh).is_some_and(|b| !b.overlaps(&opposite_box, pad)) {
            candidates.push((lh, lp));
        }
    }
    if candidates.is_empty() {
        return None;
    }
    let peeled: std::collections::HashSet<usize> = candidates
        .iter()
        .flat_map(|(_, lp)| lp.iter().copied())
        .collect();
    let rest: Vec<usize> = he_indices
        .iter()
        .copied()
        .filter(|h| !peeled.contains(h))
        .collect();
    // A peeled loop must be node-disjoint from everything the walk
    // still sees, in either direction of travel.
    let rest_nodes: std::collections::HashSet<usize> = rest
        .iter()
        .flat_map(|&h| {
            let he = &graph.half_edges[h];
            [he.origin, graph.half_edges[he.twin].origin]
        })
        .collect();
    let loops: Vec<Vec<usize>> = candidates
        .into_iter()
        .map(|(_, lp)| lp)
        .filter(|lp| {
            lp.iter()
                .all(|&h| !rest_nodes.contains(&graph.half_edges[h].origin))
        })
        .collect();
    if loops.is_empty() || rest.is_empty() {
        return None;
    }
    let peeled: std::collections::HashSet<usize> = loops.iter().flatten().copied().collect();
    let rest = he_indices
        .iter()
        .copied()
        .filter(|h| !peeled.contains(h))
        .collect();
    Some(PeeledHoles { loops, rest })
}

/// Conservative 3D box of a stored loop: its vertices plus each edge
/// curve's conservative extent.
fn loop_box(brep: &BRep, lh: LoopHandle) -> Option<Aabb> {
    let start = brep.loops.get(lh)?.half_edge;
    let mut points = Vec::new();
    let mut he = start;
    for _ in 0..=brep.half_edges.len() {
        let h = brep.half_edges.get(he)?;
        points.push(brep.vertices.get(h.origin)?.position);
        if let Some(c) = brep.edges.get(h.edge).and_then(|e| e.curve.as_ref()) {
            crate::face_box::conservative_curve_extent(c, &mut points);
        }
        he = h.next;
        if he == start {
            return Aabb::from_points(&points);
        }
    }
    None
}

/// The walk's inner-loop order key: loops of original half-edges
/// before loops carrying a record piece, then by seed (lowest index).
fn inner_key(graph: &IntersectionGraph, lp: &[usize]) -> (bool, usize) {
    let record = lp
        .iter()
        .any(|&h| matches!(graph.half_edges[h].source, EdgeSource::Shared { .. }));
    (record, lp.iter().copied().min().unwrap_or(usize::MAX))
}

/// Conservative box of a graph loop: node positions plus each
/// half-edge curve's conservative extent.
fn graph_loop_box(graph: &IntersectionGraph, lp: &[usize]) -> Option<Aabb> {
    let mut points = Vec::new();
    for &h in lp {
        let he = &graph.half_edges[h];
        points.push(graph.nodes[he.origin].position);
        let curve = he
            .curve
            .as_ref()
            .or_else(|| graph.half_edges[he.twin].curve.as_ref());
        match curve {
            Some(c) => crate::face_box::conservative_curve_extent(c, &mut points),
            // A curveless half-edge is a chord between its nodes.
            None => points.push(graph.nodes[graph.half_edges[he.twin].origin].position),
        }
    }
    Aabb::from_points(&points)
}

/// Give each peeled hole back to the reduced walk's fragment that can
/// enclose it. On a plane a fragment encloses a hole only if its outer
/// loop's box contains the hole's points, so the fragment is found by
/// box containment of the hole's first node; exactly one candidate is
/// required. The receiving fragments list their inner loops in the
/// walk's order key (the whole-face walk lists the same loops in
/// discovery order; inner-loop order carries no meaning). Returns
/// `false` (fragments untouched) when any hole has no unique fragment.
fn reattach_peeled(
    graph: &IntersectionGraph,
    fragments: &mut [Fragment],
    peel: &PeeledHoles,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> bool {
    if fragments.is_empty() {
        return false;
    }
    let boxes: Vec<_> = fragments
        .iter()
        .map(|f| graph_loop_box(graph, &f.outer_loop))
        .collect();
    let pad = 100.0 * tolerance.linear();
    let mut owner: Vec<usize> = Vec::with_capacity(peel.loops.len());
    for lp in &peel.loops {
        let Some(&first) = lp.first() else {
            return false;
        };
        let p = graph.nodes[graph.half_edges[first].origin].position;
        let point = Aabb { min: p, max: p };
        let mut candidates = boxes
            .iter()
            .enumerate()
            .filter(|(_, b)| b.as_ref().is_none_or(|b| contains(b, &point, pad)));
        match (candidates.next(), candidates.next()) {
            (Some((i, _)), None) => owner.push(i),
            _ => return false,
        }
    }
    for (lp, &i) in peel.loops.iter().zip(&owner) {
        fragments[i].inner_loops.push(lp.clone());
    }
    for f in fragments.iter_mut() {
        f.inner_loops.sort_by_key(|lp| inner_key(graph, lp));
    }
    true
}

/// `inner` lies inside `outer` grown by `pad`.
fn contains(outer: &Aabb, inner: &Aabb, pad: f64) -> bool {
    outer.min.x - pad <= inner.min.x
        && outer.min.y - pad <= inner.min.y
        && outer.min.z - pad <= inner.min.z
        && outer.max.x + pad >= inner.max.x
        && outer.max.y + pad >= inner.max.y
        && outer.max.z + pad >= inner.max.z
}

/// Compare-mode check that peeling reproduced the whole-face walk.
fn compare_peeled(peeled: &[Fragment], walked: &[Fragment]) {
    // Inner-loop ORDER is the walk's discovery order (seed order with
    // the outer swap-removed), which a merge cannot reproduce; the
    // loops themselves, their rotation and the outer must match.
    let sorted = |f: &Fragment| {
        let mut inner = f.inner_loops.clone();
        inner.sort();
        inner
    };
    let same = peeled.len() == walked.len()
        && peeled.iter().zip(walked).all(|(p, w)| {
            p.face == w.face && p.outer_loop == w.outer_loop && sorted(p) == sorted(w)
        });
    assert!(
        same,
        "peeled-hole walk differs from the whole-face walk on {:?}: peeled {:?}; walk {:?}",
        walked.first().map(|f| f.face),
        peeled
            .iter()
            .map(|f| (&f.outer_loop, &f.inner_loops))
            .collect::<Vec<_>>(),
        walked
            .iter()
            .map(|f| (&f.outer_loop, &f.inner_loops))
            .collect::<Vec<_>>()
    );
}
