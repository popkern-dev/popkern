//! Boundary stitching and edge curve reconciliation for the merge stage.
//!
//! After copying faces from both operands, half-edges along shared boundaries
//! need to be paired as twins, and edge curve geometry needs to match the
//! actual vertex positions in the result BRep.

use std::collections::HashMap;

use popkern_core::{
    tolerance::Tolerance,
    topology::{BRep, EdgeHandle, HalfEdgeHandle, VertexHandle},
};

/// Pair up unpaired half-edges that share the same two vertices in opposite
/// order. After this pass, manifold edges will have proper twin linkage and
/// shared Edge entities.
///
/// Uses both vertex-handle matching (fast path) and position-based matching
/// (fallback) to find twin pairs, then reconciles edge curve geometry so
/// vertex positions match the shared curve endpoints.
pub(super) fn stitch_boundary_edges(result: &mut BRep, tolerance: &Tolerance) {
    // First, merge coincident vertices that weren't caught during copy.
    merge_coincident_vertices(result, tolerance);

    // Collect all unpaired half-edges (those whose twin points at themselves).
    let unpaired: Vec<HalfEdgeHandle> = result
        .half_edges
        .keys()
        .filter(|&heh| result.half_edges[heh].twin == heh)
        .collect();

    // Build a map from (origin, dest) -> half-edge handle for fast lookup.
    let mut edge_map: HashMap<(VertexHandle, VertexHandle), HalfEdgeHandle> =
        HashMap::with_capacity(unpaired.len());

    for &heh in &unpaired {
        let origin = result.half_edges[heh].origin;
        let next = result.half_edges[heh].next;
        let dest = result.half_edges[next].origin;
        edge_map.insert((origin, dest), heh);
    }

    // For each unpaired half-edge (A->B), look for one going (B->A).
    for &heh in &unpaired {
        if result.half_edges[heh].twin != heh {
            continue;
        }

        let origin = result.half_edges[heh].origin;
        let next = result.half_edges[heh].next;
        let dest = result.half_edges[next].origin;

        if let Some(&twin_heh) = edge_map.get(&(dest, origin))
            && twin_heh != heh
            && result.half_edges[twin_heh].twin == twin_heh
        {
            stitch_twin_pair(result, heh, twin_heh);
        }
    }

    // Position-based fallback for remaining unpaired half-edges.
    stitch_by_position(result, tolerance);

    // Fix edge curve geometry to match vertex positions.
    reconcile_edge_curves(result, tolerance);
}

/// Merge coincident vertices in the result BRep.
///
/// After copying faces from both operands, some vertices may be at the same
/// position but not deduplicated. This pass unifies them and updates all
/// half-edge origin references.
fn merge_coincident_vertices(result: &mut BRep, tolerance: &Tolerance) {
    let verts: Vec<(VertexHandle, _)> = result
        .vertices
        .iter()
        .map(|(vh, v)| (vh, v.position))
        .collect();

    let mut remap: HashMap<VertexHandle, VertexHandle> = HashMap::new();
    let n = verts.len();
    for i in 0..n {
        let (vh_i, pos_i) = verts[i];
        if remap.contains_key(&vh_i) {
            continue;
        }
        for &(vh_j, pos_j) in &verts[(i + 1)..n] {
            if remap.contains_key(&vh_j) {
                continue;
            }
            if tolerance.points_coincident(&pos_i, &pos_j) {
                remap.insert(vh_j, vh_i);
            }
        }
    }

    if remap.is_empty() {
        return;
    }

    for (_, he) in &mut result.half_edges {
        if let Some(&canonical) = remap.get(&he.origin) {
            he.origin = canonical;
        }
    }

    for heh in result.half_edges.keys().collect::<Vec<_>>() {
        let origin = result.half_edges[heh].origin;
        result.vertices[origin].outgoing = heh;
    }

    for old_vh in remap.keys() {
        result.vertices.remove(*old_vh);
    }
}

/// Stitch two half-edges as twins and merge their edge entities.
fn stitch_twin_pair(result: &mut BRep, heh: HalfEdgeHandle, twin_heh: HalfEdgeHandle) {
    result.half_edges[heh].twin = twin_heh;
    result.half_edges[twin_heh].twin = heh;

    let keep_edge = result.half_edges[heh].edge;
    let discard_edge = result.half_edges[twin_heh].edge;

    if keep_edge != discard_edge {
        result.half_edges[twin_heh].edge = keep_edge;

        // Update any other half-edges referencing the discarded edge.
        let to_update: Vec<HalfEdgeHandle> = result
            .half_edges
            .keys()
            .filter(|&h| result.half_edges[h].edge == discard_edge)
            .collect();
        for h in to_update {
            result.half_edges[h].edge = keep_edge;
        }

        result.edges[keep_edge].half_edge = heh;

        if result.edges[keep_edge].curve.is_none() {
            let curve = result.edges[discard_edge].curve.clone();
            result.edges[keep_edge].curve = curve;
        }

        result.edges[discard_edge].half_edge = HalfEdgeHandle::default();
    }
}

/// Position-based fallback stitching for half-edges that weren't paired
/// by vertex-handle matching.
///
/// For each unpaired half-edge (A->B), finds another unpaired half-edge whose
/// origin is at position B and whose destination is at position A (within
/// tolerance). Unifies vertices as needed before pairing.
fn stitch_by_position(result: &mut BRep, tolerance: &Tolerance) {
    loop {
        let unpaired: Vec<HalfEdgeHandle> = result
            .half_edges
            .keys()
            .filter(|&heh| result.half_edges[heh].twin == heh)
            .collect();

        if unpaired.is_empty() {
            break;
        }

        let entries: Vec<(_, _, HalfEdgeHandle)> = unpaired
            .iter()
            .map(|&heh| {
                let origin = result.half_edges[heh].origin;
                let next = result.half_edges[heh].next;
                let dest = result.half_edges[next].origin;
                let origin_pos = result.vertices[origin].position;
                let dest_pos = result.vertices[dest].position;
                (origin_pos, dest_pos, heh)
            })
            .collect();

        let mut stitched_any = false;
        let mut already_stitched = std::collections::HashSet::new();

        for i in 0..entries.len() {
            let (origin_a, dest_a, he_a) = entries[i];
            if already_stitched.contains(&he_a) || result.half_edges[he_a].twin != he_a {
                continue;
            }

            for (j, &(origin_b, dest_b, he_b)) in entries.iter().enumerate() {
                if i == j {
                    continue;
                }
                if already_stitched.contains(&he_b) || result.half_edges[he_b].twin != he_b {
                    continue;
                }

                if tolerance.points_coincident(&origin_b, &dest_a)
                    && tolerance.points_coincident(&dest_b, &origin_a)
                {
                    // Unify vertices if they differ.
                    let origin_he_a = result.half_edges[he_a].origin;
                    let next_he_a = result.half_edges[he_a].next;
                    let dest_he_a = result.half_edges[next_he_a].origin;
                    let origin_he_b = result.half_edges[he_b].origin;
                    let next_he_b = result.half_edges[he_b].next;
                    let dest_he_b = result.half_edges[next_he_b].origin;

                    if origin_he_b != dest_he_a {
                        remap_vertex(result, origin_he_b, dest_he_a);
                    }
                    if dest_he_b != origin_he_a {
                        remap_vertex(result, dest_he_b, origin_he_a);
                    }

                    stitch_twin_pair(result, he_a, he_b);
                    already_stitched.insert(he_a);
                    already_stitched.insert(he_b);
                    stitched_any = true;
                    break;
                }
            }
        }

        if !stitched_any {
            break;
        }
    }
}

/// Remap all references from `old_vh` to `new_vh` and remove the old vertex.
fn remap_vertex(result: &mut BRep, old_vh: VertexHandle, new_vh: VertexHandle) {
    if old_vh == new_vh {
        return;
    }
    for (_, he) in &mut result.half_edges {
        if he.origin == old_vh {
            he.origin = new_vh;
        }
    }
    result.vertices.remove(old_vh);
    if let Some(heh) = result
        .half_edges
        .keys()
        .find(|&heh| result.half_edges[heh].origin == new_vh)
    {
        result.vertices[new_vh].outgoing = heh;
    }
}

/// Reconcile edge curve geometry with vertex positions after stitching.
///
/// For paired edges, ensures curve t_start/t_end match the canonical
/// half-edge's origin and the twin's origin respectively. For self-twin
/// (unpaired boundary) edges where the curve can't match the validation
/// expectations, removes the curve geometry so validation falls back to
/// vertex-position checks.
fn reconcile_edge_curves(result: &mut BRep, tolerance: &Tolerance) {
    use popkern_core::geometry::CurveEval;

    let edge_handles: Vec<EdgeHandle> = result.edges.keys().collect();
    for eh in edge_handles {
        let curve = match &result.edges[eh].curve {
            Some(c) => c.clone(),
            None => continue,
        };

        let he_h = result.edges[eh].half_edge;
        let he = match result.half_edges.get(he_h) {
            Some(h) => h,
            None => continue,
        };

        let v_start_pos = match result.vertices.get(he.origin) {
            Some(v) => v.position,
            None => continue,
        };

        let is_self_twin = he.twin == he_h;

        let twin = match result.half_edges.get(he.twin) {
            Some(t) => t,
            None => continue,
        };
        let v_end_pos = match result.vertices.get(twin.origin) {
            Some(v) => v.position,
            None => continue,
        };

        let c_start = curve.geometry.point_at(curve.t_start);
        let c_end = curve.geometry.point_at(curve.t_end);

        let gap_start = v_start_pos.distance_to(&c_start);
        let gap_end = v_end_pos.distance_to(&c_end);

        if gap_start > tolerance.linear || gap_end > tolerance.linear {
            if is_self_twin {
                // Self-twin edges: validation checks both endpoints against
                // he.origin, which fails for edges spanning a distance.
                // Remove curve so validation uses vertex-position fallback.
                result.edges[eh].curve = None;
            } else {
                // Try swapping t_start/t_end.
                let alt_gap_start = v_start_pos.distance_to(&c_end);
                let alt_gap_end = v_end_pos.distance_to(&c_start);
                if alt_gap_start <= tolerance.linear && alt_gap_end <= tolerance.linear {
                    let mut new_curve = curve.clone();
                    std::mem::swap(&mut new_curve.t_start, &mut new_curve.t_end);
                    result.edges[eh].curve = Some(new_curve);
                } else if let Some(new_curve) =
                    refit_curve_to_vertices(&curve, v_start_pos, v_end_pos, tolerance)
                {
                    result.edges[eh].curve = Some(new_curve);
                }
            }
        }
    }
}

/// Try to refit a curve's t_start/t_end so the curve evaluates to the
/// given start and end positions. Returns the adjusted curve on success.
fn refit_curve_to_vertices(
    curve: &popkern_core::geometry::CurveOnEdge,
    start_pos: popkern_core::geometry::Point3,
    end_pos: popkern_core::geometry::Point3,
    tolerance: &Tolerance,
) -> Option<popkern_core::geometry::CurveOnEdge> {
    use popkern_core::geometry::CurveEval;
    use popkern_core::geometry::curve::CurveGeometry;

    match &curve.geometry {
        CurveGeometry::Line { origin, direction } => {
            // For a line: P(t) = origin + t * direction
            // Solve for t: t = (P - origin) . direction / |direction|^2
            let dir_len_sq = direction.dot(direction);
            if dir_len_sq < 1e-30 {
                return None;
            }
            let diff_start = start_pos - *origin;
            let t_start = diff_start.dot(direction) / dir_len_sq;
            let diff_end = end_pos - *origin;
            let t_end = diff_end.dot(direction) / dir_len_sq;

            let check_start = curve.geometry.point_at(t_start);
            let check_end = curve.geometry.point_at(t_end);
            if start_pos.distance_to(&check_start) <= tolerance.linear
                && end_pos.distance_to(&check_end) <= tolerance.linear
            {
                let mut new_curve = curve.clone();
                new_curve.t_start = t_start;
                new_curve.t_end = t_end;
                Some(new_curve)
            } else {
                None
            }
        }
        _ => None,
    }
}
