//! Cellular boolean: the primary boolean pipeline using cell decomposition.
//!
//! Reuses the intersection graph infrastructure (phases 1-5 of graph-merge)
//! to produce face fragments, then classifies fragments via point-in-solid
//! to determine cell membership, selects cells by operation predicate, and
//! derives fragment selections from the cell selection.
//!
//! The assembly step reuses graph-merge's `assemble_result`, so the only
//! difference from graph-merge is HOW fragments are selected — cellular
//! selection is order-independent because it classifies volumetric cells
//! rather than individual fragments.

use super::record_rotation::rotate_closed_record_to_hits;

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::complex::serve::SelectionSource;
use crate::error::{BoolDiagnostics, BoolError, BoolOperation};
use crate::intersect::IntersectionResult;
use crate::merge::graph::{IntersectionGraph, Operand};
use crate::merge::graph_duplicates::remove_collinear_duplicates;
use crate::merge::graph_intersect::insert_intersection_curves_at;
use crate::merge::graph_split::split_shared_edges_at_nodes;
use crate::merge::graph_split_shared::split_edges_at_shared_nodes;
use crate::merge::radial::sort_all_radial;
use crate::merge::walk::walk_all_fragments_ext;

/// Execute a boolean operation using the cellular decomposition approach.
///
/// Phases 1-5 (graph build, edge-face splits, intersection curves, radial
/// sort, fragment walk) are identical to graph-merge. Phase 6 replaces
/// per-fragment ray casting with cell-based classification:
///
/// 1. For each fragment, classify a sample point against both bodies
/// 2. Derive a cell membership (in A, in B, in both, in neither)
/// 3. Select cells by the operation predicate
/// 4. Keep fragments whose cell is selected
///
/// This is order-independent because the classification depends only on
/// geometry (point-in-solid), not on iteration order.
///
/// Returns the body and its CARRIED edges — those built from unsplit
/// operand edges, whose curve-surface agreement the operand already
/// had (the post-check measures the rest).
pub(crate) fn cellular_graph_merge(
    original_a: &BRep,
    original_b: &BRep,
    intersection_result: &IntersectionResult,
    op: BoolOperation,
    tolerance: &Tolerance,
    scopes: Option<&crate::scope::OperandScopes>,
    source: SelectionSource,
) -> Result<(BRep, crate::serve_region::Carried), BoolDiagnostics> {
    let snap_tol = tolerance.linear() * 2.0;
    let mut graph = IntersectionGraph::new();

    // Process edge-face hits.
    if popkern_core::trace::trace_enabled("boolean") {
        eprintln!(
            "  [cellular] {} edge-face hits",
            intersection_result.edge_face_hits.len()
        );
        for hit in &intersection_result.edge_face_hits {
            eprintln!(
                "    hit at ({:.3},{:.3},{:.3})",
                hit.point.x, hit.point.y, hit.point.z
            );
        }
    }
    let operand_pair = [original_a, original_b];
    let breps = Some((&operand_pair[..], tolerance.linear()));

    // Insert intersection curves as shared edges, then
    // split them at edge-face hit nodes they pass through (the loop
    // walk must be able to route through those junctions).
    // Closed intersection curves must CLOSE at an edge-face hit
    // point when one lies on them: the loop walk enters/bounces off
    // an intersection circle only at its closure node (§7.4), and
    // the analytic cylinder-hole cases have always closed at the
    // cylinder seam — the crossing — by construction. A circle from
    // an offset walk start (or a swept boss whose seam is elsewhere)
    // must be rotated to that convention.
    // With the table deciding vertices, it decides where closed
    // records close (`vertex_identity::decide`).
    let table_enabled = crate::tolerance_zone::table::table_enabled();
    let mut rotated: Vec<crate::intersect::IntersectionRecord> = intersection_result
        .records
        .iter()
        .map(|r| match table_enabled {
            true => r.clone(),
            false => {
                rotate_closed_record_to_hits(r, &intersection_result.edge_face_hits, tolerance)
            }
        })
        .collect();
    // COINCIDENT-INTERFACE containment (PK-465): trim residues are
    // unproven on booleans that carry coplanar face overlaps — the
    // tuned flush arbitration owns those interfaces, and residues
    // riding them destabilized the flush-gusset family's walks
    // (29.6 loop gap measured). Strip every residue there; the
    // transversal fixtures the mechanism serves (the angled tap)
    // carry no overlaps.
    if !intersection_result.coincidence.face_overlaps.is_empty() {
        rotated.retain(|r| r.residue_owner_a.is_none());
    }
    crate::intersect_residue_filter::drop_kept_flank_residues(
        &mut rotated,
        op,
        original_a,
        original_b,
        &intersection_result.aabbs_a,
        &intersection_result.aabbs_b,
        tolerance,
    );
    // Multi-turn bands get a node where a record crosses their period
    // lattice above the seam (PK-919).
    let rotated = crate::intersect_band_lattice::split_records_at_band_lattice(
        rotated, original_a, original_b,
    );
    let partial_overlaps = intersection_result
        .coincidence
        .face_overlaps
        .iter()
        .any(|o| !o.full_overlap);
    let (rotated, identity) = super::vertex_identity::decide(
        (original_a, original_b),
        intersection_result,
        rotated,
        tolerance,
        partial_overlaps,
    );
    let graph_splits_junctions = identity.graph_identity();

    // Insert original topology from both operands, each table-held
    // vertex at its anchor.
    // A scoped run (PK-814, PK-816) builds each scoped operand's graph
    // from its active faces only; the assembly carries the far faces the
    // operation keeps and the rest never form fragments.
    let [anchors_a, anchors_b] = &identity.vertex_anchors;
    let scope = scopes.and_then(|s| s.a.as_ref());
    let (_vmap_a, hemap_a) = crate::merge::graph_build::insert_operand_topology_scoped(
        &mut graph,
        original_a,
        Operand::A,
        snap_tol,
        scope,
        Some(anchors_a),
    );
    let (_vmap_b, hemap_b) = crate::merge::graph_build::insert_operand_topology_scoped(
        &mut graph,
        original_b,
        Operand::B,
        snap_tol,
        scopes.and_then(|s| s.b.as_ref()),
        Some(anchors_b),
    );
    crate::merge::graph_split_node::split_edges_at_face_hits(
        &mut graph,
        &identity.hits,
        snap_tol,
        breps,
    );
    insert_intersection_curves_at(&mut graph, &rotated, snap_tol, identity.ends.as_deref());
    // Split SHARED edges at record×record T-junctions: one record's
    // endpoint frequently lands mid-span of another (a dovetail
    // rail's flank trace ends ON the joint-shoulder trace), and the
    // walk cannot route the junction without a node-split there
    // (measured: the rail's floor strip never became a fragment and
    // the carve went non-manifold). Every node is a candidate — the
    // splitter skips nodes not on a shared edge's interior. With the
    // table deciding vertices the records arrived split at its paves
    // (T-junction consequences, plan §4.3a), common-block scenes too,
    // and the pass serves only the graph's own identity: a declined
    // table or the kill switch.
    if graph_splits_junctions {
        let all_nodes: Vec<usize> = (0..graph.nodes.len()).collect();
        let before = graph.half_edges.len();
        split_shared_edges_at_nodes(&mut graph, &all_nodes, snap_tol);
        crate::complex_shadow::subsumed_pass_fired(
            "split_shared_edges_at_nodes",
            (graph.half_edges.len() - before) / 2,
            0.0,
        );
    }
    if popkern_core::trace::trace_enabled("boolean") {
        let mut i = 0;
        while i < graph.half_edges.len() {
            let he = &graph.half_edges[i];
            if i < he.twin && matches!(he.source, crate::merge::graph::EdgeSource::Shared { .. }) {
                let a = graph.nodes[he.origin].position;
                let b = graph.nodes[graph.half_edges[he.twin].origin].position;
                eprintln!(
                    "  [graph] shared he {i}: ({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3}) curve={}",
                    a.x,
                    a.y,
                    a.z,
                    b.x,
                    b.y,
                    b.z,
                    he.curve.is_some()
                );
            }
            i += 2;
        }
    }

    // Phase 3.1–3.2: Handle coincident face boundaries; also split
    // graph edges at nodes created by edge-face hits — an
    // intersection curve passing THROUGH such a node (e.g. a seam
    // edge crossing a boss/plate intersection circle) must get a
    // node there or the loop walk cannot route through it.
    if !intersection_result.edge_face_hits.is_empty() || partial_overlaps {
        // With the table deciding, the records' ends split the operand
        // edges it put them on; the dead-end completion (march slop
        // beyond the table's reach) still runs.
        if graph_splits_junctions {
            let before = graph.half_edges.len();
            split_edges_at_shared_nodes(&mut graph, snap_tol, breps);
            crate::complex_shadow::subsumed_pass_fired(
                "split_edges_at_shared_nodes",
                graph.half_edges.len().saturating_sub(before) / 2,
                0.0,
            );
        } else {
            crate::merge::graph_paves::split_operand_edges_at_paves(
                &mut graph,
                &identity.operand_paves,
                snap_tol,
                breps,
            );
            crate::merge::graph_split_shared::complete_dead_end_cuts(&mut graph, snap_tol, breps);
        }
        // A record along an operand edge is one edge with it: the
        // table's common blocks decide which (plan §4.3a). The graph's
        // collinear-duplicate search follows, deciding alone under its
        // own identity and, on the table's path, the residue: curves
        // the table keeps apart (beyond its vertices' tolerances) that
        // the walk's own uv dedup still merges on one face, within
        // 100 × linear (measured: a near-tangent cap's 1e-5 lens).
        if partial_overlaps && !graph_splits_junctions {
            crate::merge::graph_duplicates::remove_common_block_duplicates(
                &mut graph,
                &identity.common_blocks,
                &rotated,
                snap_tol,
            );
        }
        if partial_overlaps && graph_splits_junctions {
            let before = graph.duplicate_of.len();
            remove_collinear_duplicates(&mut graph, tolerance.linear());
            crate::complex_shadow::subsumed_pass_fired(
                "remove_collinear_duplicates",
                graph.duplicate_of.len() - before,
                0.0,
            );
        }
    }

    if crate::complex_shadow::shadow_enabled() {
        crate::tolerance_zone::shadow::run(
            &graph,
            (original_a, original_b),
            intersection_result,
            &rotated,
            tolerance,
        );
    }

    // Sort radial edges at all nodes.
    sort_all_radial(&mut graph);

    // Walk loops to find face fragments.
    let overlap_pairs: Vec<_> = intersection_result
        .coincidence
        .face_overlaps
        .iter()
        .map(|o| {
            let face = |operand, face| crate::merge::graph::FaceRef { operand, face };
            (face(Operand::A, o.face_a), face(Operand::B, o.face_b))
        })
        .collect();
    let mut walk_fallbacks: Vec<(crate::merge::graph::FaceRef, String)> = Vec::new();
    let holes_of = |s: Option<&crate::scope::ActiveScope>| -> Option<std::collections::HashSet<_>> {
        s.map(|s| s.far_holes.values().flatten().copied().collect())
    };
    let carried_holes = [holes_of(scope), holes_of(scopes.and_then(|s| s.b.as_ref()))];
    let fragments = walk_all_fragments_ext(
        &graph,
        &[original_a, original_b],
        !intersection_result.coincidence.face_overlaps.is_empty(),
        true,
        &overlap_pairs,
        &mut walk_fallbacks,
        Some(&crate::merge::walk::OperandHalfEdgeMaps {
            a: &hemap_a,
            b: &hemap_b,
            boxes_a: &intersection_result.aabbs_a,
            boxes_b: &intersection_result.aabbs_b,
            carried_holes: [carried_holes[0].as_ref(), carried_holes[1].as_ref()],
        }),
        tolerance,
    );

    for (face_ref, reason) in &walk_fallbacks {
        crate::chart_walk_fallback::record_chart_walk_fallback(face_ref, reason);
    }
    // A face's own hole the walk filed apart rejoins the piece around
    // it, for every classifier: served as a piece it is an inverted
    // patch over a face that lost its hole.
    let fragments = crate::complex::hole_phantoms::rejoin_own_holes(
        &graph,
        &fragments,
        &[original_a, original_b],
        tolerance.linear(),
    )
    .unwrap_or(fragments);

    let frontier = frontier_keys(
        &graph,
        [original_a, original_b],
        [scope, scopes.and_then(|s| s.b.as_ref())],
    );
    let complex_input = crate::complex::pipeline::ComplexInput {
        graph: &graph,
        fragments: &fragments,
        operands: &[original_a, original_b],
        boxes: &[
            &intersection_result.aabbs_a[..],
            &intersection_result.aabbs_b[..],
        ],
        overlaps: &overlap_pairs,
        tolerance,
        frontier: &frontier,
    };
    let cellular = || -> Result<_, BoolDiagnostics> {
        // Phase 6 (CELLULAR): Classify fragments by cell membership.
        let mut selections = super::classify_fragments::classify_fragments_cellular(
            &fragments,
            &graph,
            original_a,
            original_b,
            op,
            tolerance,
            &intersection_result.coincidence.face_overlaps,
            (&intersection_result.aabbs_a, &intersection_result.aabbs_b),
        )?;
        // A mixed-sense pair wires as a same-sense pair (PK-913): B's
        // kept fragments flip when every shared record says it winds
        // opposite to A, and its carried far faces reverse with them
        // (PK-916).
        let (winding, flipped_b) =
            super::operand_winding::agree_operand_winding(&graph, &fragments, &mut selections);
        if crate::complex_shadow::shadow_enabled() {
            crate::complex::shadow_compare::shadow_compare(&complex_input, &selections, op);
        }
        Ok((selections, winding, flipped_b))
    };
    // The cell complex decides every fragment (PK-976): region
    // membership by propagation, one skin per same-domain face, one
    // folded-normal convention (no winding vote).
    let mut split = None;
    let (selections, winding, flipped_b) = match source {
        SelectionSource::Complex | SelectionSource::ComplexSplit => {
            let (selections, plan) = crate::complex::serve::complex_selections(
                &complex_input,
                op,
                source == SelectionSource::ComplexSplit,
            )?;
            split = plan;
            let flipped_b = crate::complex::serve::b_reversed(&fragments, &selections, op);
            (selections, Default::default(), flipped_b)
        }
        SelectionSource::ComplexFirst => {
            match crate::complex::serve::complex_selections(&complex_input, op, false) {
                Ok((selections, _)) => {
                    let flipped_b = crate::complex::serve::b_reversed(&fragments, &selections, op);
                    (selections, Default::default(), flipped_b)
                }
                Err(diag) if matches!(diag.error, Some(BoolError::NonManifoldResult { .. })) => {
                    return Err(diag);
                }
                Err(diag) => {
                    crate::complex_shadow::record_complex_primary_fallback(op, &diag);
                    cellular()?
                }
            }
        }
        SelectionSource::Cellular => cellular()?,
    };
    let (far_a, far_b) = scopes.map(|s| s.carried(op)).unwrap_or_default();
    if popkern_core::trace::trace_enabled("boolean") && winding.same > 0 {
        eprintln!(
            "  [cellular] shared winding: {} same, {} opposite; B {}",
            winding.same,
            winding.opposite,
            if flipped_b && far_b.is_some() {
                "flipped to agree with A, its carried far faces reversed"
            } else if flipped_b {
                "flipped to agree with A"
            } else {
                "left"
            }
        );
    }

    // Order-independent assembly from selected fragments.
    let kept: Vec<usize> = (0..fragments.len())
        .filter(|&i| selections[i].keep)
        .collect();

    let mut carried = crate::serve_region::Carried::default();
    let result = super::assemble::cellular_assemble(
        &graph,
        &fragments,
        &selections,
        &[original_a, original_b],
        *tolerance,
        &mut carried,
        super::assemble_far_carry::FarCarry {
            a: far_a,
            b: far_b,
            b_reversed: flipped_b,
        },
        split.as_ref(),
    )
    .map_err(|mut diag| {
        // ChartSubdivisionReport promotion (PK-106): a failed cellular
        // boolean carries, per face that could not chart-walk, WHY it
        // took the 3D-Newell fallback — the substrate defects are usually
        // the trailhead for the assembly failure being reported.
        for (face_ref, reason) in &walk_fallbacks {
            diag.add(crate::error::DiagnosticEntry {
                severity: crate::error::Severity::Warning,
                phase: crate::error::Stage::Merge,
                message: format!(
                    "chart walk fell back to the 3D-Newell walk on {face_ref:?}: {reason}"
                ),
                entities: vec![crate::error::EntityRef::Face(face_ref.face)],
                measured_value: None,
                location: None,
                remediation: None,
            });
        }
        diag
    })?;

    // Build provenance map (available via boolean_with_provenance).
    let _provenance =
        super::assemble_cleanup::build_face_provenance(&result, &kept, &fragments, &selections);

    Ok((result, carried))
}

/// Graph edge keys (`min(half-edge, twin)`) of the edges bounding
/// each operand's scope, split pieces included: one side's half-edge
/// entered the graph, the other's did not, so the complex sees one fin
/// there by design. That is an active face against a far face, or
/// against a far hole of an active face (measured: a pipped die's
/// dimple rims, carried as far holes of the face they cut, around
/// dimple faces that stayed active). Read from each graph edge's
/// provenance, which names every piece of a split edge.
fn frontier_keys(
    graph: &IntersectionGraph,
    operands: [&BRep; 2],
    scopes: [Option<&crate::scope::ActiveScope>; 2],
) -> std::collections::BTreeSet<usize> {
    use crate::merge::graph::EdgeSource;
    let bounds = |k: usize, e: &popkern_core::topology::EdgeHandle| {
        scopes[k].is_some_and(|s| s.bounds(operands[k], *e))
    };
    graph
        .half_edges
        .iter()
        .enumerate()
        .filter(|(_, he)| match &he.source {
            EdgeSource::Original(o, e)
            | EdgeSource::Split {
                operand: o,
                original: e,
            } => bounds(o.index(), e),
            EdgeSource::Shared { .. } => false,
        })
        .map(|(h, he)| h.min(he.twin))
        .collect()
}
