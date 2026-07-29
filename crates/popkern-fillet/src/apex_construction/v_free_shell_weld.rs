//! Seam weld + vertex unification for the V-free consolidation
//! (PK-153): dead boundary sides with coincident antiparallel spans
//! re-twin pairwise (duplicate edge records unify; curve authority
//! prefers the non-linear record), and vertex records duplicated
//! across welded seams merge by position.

use std::collections::{BTreeMap, BTreeSet};

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::topology::{BRep, EdgeHandle, HalfEdgeHandle, VertexHandle};

use super::side_span;

/// Re-twin dead boundary sides whose spans coincide antiparallel,
/// unifying the duplicate edge records. Record survivorship:
/// `preferred` handles win (downstream plans hold them — retiring a
/// plan-referenced record strands every later consumer, measured as
/// the exact-boundary `stale_edge` decline under the PK-167 in-place
/// rewrite); otherwise the side carrying a non-linear curve wins.
/// When a preferred linear record beats a curved duplicate, the
/// survivor ADOPTS the curve so authority is never lost.
pub(in crate::apex_construction) fn weld_dead_coincident_sides(
    brep: &mut BRep,
    preferred: &BTreeSet<EdgeHandle>,
) -> (usize, Vec<VertexHandle>) {
    let tol = brep.tolerance.linear().max(1.0e-9) * 10.0;
    let dead: Vec<HalfEdgeHandle> = brep
        .half_edges
        .iter()
        .filter(|(_, data)| !brep.half_edges.contains_key(data.twin))
        .map(|(handle, _)| handle)
        .collect();
    let mut welded = 0usize;
    let mut seam_vertices: Vec<VertexHandle> = Vec::new();
    let mut used: std::collections::BTreeSet<HalfEdgeHandle> = std::collections::BTreeSet::new();
    for i in 0..dead.len() {
        if used.contains(&dead[i]) {
            continue;
        }
        let Some((a_start, a_end)) = side_span(brep, dead[i]) else {
            continue;
        };
        for j in (i + 1)..dead.len() {
            if used.contains(&dead[j]) {
                continue;
            }
            let Some((b_start, b_end)) = side_span(brep, dead[j]) else {
                continue;
            };
            if a_start.distance_to(&b_end) > tol || a_end.distance_to(&b_start) > tol {
                continue;
            }
            // Endpoint coincidence is not curve identity: a fillet
            // arc and its chord share endpoints but differ by the
            // sagitta at midpoint (measured 2.93e-2 for an r=0.10
            // quarter arc — welding that pair put the arc's 3D
            // curve on an edge whose plane-side pcurves draw the
            // chord, a body-level lie STEP import flags at exactly
            // that number). Both sides must also agree at their
            // curve midpoints before pairing; a curve-less side
            // contributes its chord midpoint.
            let midpoint_of = |side: HalfEdgeHandle, start: Point3, end: Point3| -> Point3 {
                brep.half_edges
                    .get(side)
                    .and_then(|he| brep.edges.get(he.edge))
                    .and_then(|edge| edge.curve.as_ref())
                    .map(|coe| {
                        use popkern_core::geometry::curve::CurveEval;
                        coe.geometry.point_at(0.5 * (coe.t_start + coe.t_end))
                    })
                    .unwrap_or_else(|| {
                        Point3::new(
                            0.5 * (start.x + end.x),
                            0.5 * (start.y + end.y),
                            0.5 * (start.z + end.z),
                        )
                    })
            };
            if midpoint_of(dead[i], a_start, a_end)
                .distance_to(&midpoint_of(dead[j], b_start, b_end))
                > tol
            {
                continue;
            }
            let (side_a, side_b) = (dead[i], dead[j]);
            let edge_a = brep.half_edges[side_a].edge;
            let edge_b = brep.half_edges[side_b].edge;
            // Curve authority: keep the non-linear (blend arc) record.
            let is_curved = |edge: popkern_core::topology::EdgeHandle| {
                brep.edges
                    .get(edge)
                    .and_then(|e| e.curve.as_ref())
                    .is_some_and(|c| !matches!(c.geometry, CurveGeometry::Line { .. }))
            };
            let b_wins = (preferred.contains(&edge_b) && !preferred.contains(&edge_a))
                || (!preferred.contains(&edge_a) && is_curved(edge_b) && !is_curved(edge_a));
            let (keep, drop, drop_side) = if edge_a == edge_b {
                (edge_a, None, side_b)
            } else if b_wins {
                (edge_b, Some(edge_a), side_a)
            } else {
                (edge_a, Some(edge_b), side_b)
            };
            // Curve-authority transplant: a preferred survivor that
            // is linear adopts the dropped record's non-linear curve
            // (and its recorded tolerance) before the duplicate dies.
            if let Some(dropped) = drop
                && is_curved(dropped)
                && !is_curved(keep)
            {
                let adopted = brep
                    .edges
                    .get(dropped)
                    .map(|e| (e.curve.clone(), e.tolerance));
                if let (Some((curve, tolerance)), Some(keep_data)) =
                    (adopted, brep.edges.get_mut(keep))
                {
                    keep_data.curve = curve;
                    if tolerance.is_some() {
                        keep_data.tolerance = tolerance;
                    }
                }
            }
            brep.half_edges[side_a].twin = side_b;
            brep.half_edges[side_b].twin = side_a;
            brep.half_edges[drop_side].edge = keep;
            brep.half_edges[drop_side].pcurve = None;
            if let Some(drop_edge) = drop {
                brep.edges.remove(drop_edge);
                brep.boundary_edges.remove(&drop_edge);
            }
            brep.boundary_edges.remove(&keep);
            if let Some(edge_data) = brep.edges.get_mut(keep) {
                edge_data.half_edge = side_a;
            }
            for side in [side_a, side_b] {
                if let Some(data) = brep.half_edges.get(side) {
                    seam_vertices.push(data.origin);
                    if let Some(next) = brep.half_edges.get(data.next) {
                        seam_vertices.push(next.origin);
                    }
                }
            }
            used.insert(side_a);
            used.insert(side_b);
            welded += 1;
            break;
        }
    }
    if crate::direct_apex_trace::enabled() {
        for side in &dead {
            if used.contains(side) || !brep.half_edges.contains_key(*side) {
                continue;
            }
            if let Some((start, end)) = side_span(brep, *side) {
                let face = brep
                    .loops
                    .get(brep.half_edges[*side].loop_)
                    .map(|l| format!("{:?}", l.face))
                    .unwrap_or_else(|| "?".into());
                eprintln!(
                    "[direct-apex-trace] consolidation: unwelded dead side {side:?} \
                     face={face} span=({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3})",
                    start.x, start.y, start.z, end.x, end.y, end.z
                );
            }
        }
    }
    (welded, seam_vertices)
}

/// Merge coincident vertex records referenced across welded seams.
/// Scoped to the weld-adjacent vertices — position-coincident distinct
/// vertices elsewhere (periodic seams) are by-design and untouched.
/// `preferred` handles win canonicalship (downstream plans hold them
/// — removing a plan-referenced record strands later consumers on an
/// invalid SlotMap key, measured at the pyramid3 exact-trim author).
pub(in crate::apex_construction) fn merge_coincident_seam_vertices(
    brep: &mut BRep,
    seam: &[VertexHandle],
    preferred: &BTreeSet<VertexHandle>,
) -> BTreeMap<VertexHandle, VertexHandle> {
    let tol = brep.tolerance.linear().max(1.0e-9) * 10.0;
    let mut unique: Vec<VertexHandle> = seam.to_vec();
    unique.sort();
    unique.dedup();
    let vertices: Vec<(VertexHandle, Point3)> = unique
        .into_iter()
        .filter_map(|handle| {
            brep.vertices
                .get(handle)
                .map(|data| (handle, data.position))
        })
        .collect();
    let mut canonical: BTreeMap<VertexHandle, VertexHandle> = BTreeMap::new();
    for i in 0..vertices.len() {
        for j in (i + 1)..vertices.len() {
            if vertices[i].1.distance_to(&vertices[j].1) <= tol {
                let target = *canonical.get(&vertices[i].0).unwrap_or(&vertices[i].0);
                canonical.entry(vertices[j].0).or_insert(target);
            }
        }
    }
    // Preferred canonicalship: re-root every cluster whose members
    // include a preferred handle so the plan-referenced record
    // survives the merge.
    if !preferred.is_empty() {
        let mut cluster_root: BTreeMap<VertexHandle, VertexHandle> = BTreeMap::new();
        for (&dup, &target) in &canonical {
            let root = *cluster_root.get(&target).unwrap_or(&target);
            if preferred.contains(&dup) && !preferred.contains(&root) {
                cluster_root.insert(target, dup);
            }
        }
        if !cluster_root.is_empty() {
            let mut rerooted: BTreeMap<VertexHandle, VertexHandle> = BTreeMap::new();
            for (&dup, &target) in &canonical {
                let root = *cluster_root.get(&target).unwrap_or(&target);
                if dup == root {
                    // The new root must not be removed; its old root
                    // becomes the duplicate instead.
                    if target != root {
                        rerooted.insert(target, root);
                    }
                } else {
                    rerooted.insert(dup, root);
                }
            }
            canonical = rerooted;
        }
    }
    if canonical.is_empty() {
        return BTreeMap::new();
    }
    if std::env::var_os("PK_SLIR_TRACE").is_some() {
        for (dup, target) in &canonical {
            eprintln!(
                "[slir] vertex merge {dup:?} -> {target:?} (dup preferred: {}, target preferred: {})",
                preferred.contains(dup),
                preferred.contains(target)
            );
        }
    }
    let handles: Vec<HalfEdgeHandle> = brep.half_edges.keys().collect();
    for handle in handles {
        let origin = brep.half_edges[handle].origin;
        if let Some(target) = canonical.get(&origin) {
            brep.half_edges[handle].origin = *target;
        }
    }
    let mut merged: BTreeMap<VertexHandle, VertexHandle> = BTreeMap::new();
    for (duplicate, target) in canonical {
        if brep.vertices.remove(duplicate).is_some() {
            merged.insert(duplicate, target);
        }
    }
    merged
}
