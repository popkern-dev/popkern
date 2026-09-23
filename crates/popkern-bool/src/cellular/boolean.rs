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

use crate::error::{BoolDiagnostics, BoolOperation};
use crate::intersect::IntersectionResult;
use crate::merge::graph::{IntersectionGraph, Operand};
use crate::merge::graph_build::insert_operand_topology;
use crate::merge::graph_intersect::insert_intersection_curves;
use crate::merge::graph_split::{remove_collinear_duplicates, split_shared_edges_at_nodes};
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
pub(crate) fn cellular_graph_merge(
    original_a: &BRep,
    original_b: &BRep,
    intersection_result: &IntersectionResult,
    op: BoolOperation,
    tolerance: &Tolerance,
) -> Result<BRep, BoolDiagnostics> {
    let snap_tol = tolerance.linear() * 2.0;
    let mut graph = IntersectionGraph::new();

    // Insert original topology from both operands.
    let (_vmap_a, _hemap_a) = insert_operand_topology(&mut graph, original_a, Operand::A, snap_tol);
    let (_vmap_b, _hemap_b) = insert_operand_topology(&mut graph, original_b, Operand::B, snap_tol);

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
    let breps = Some((original_a, original_b, tolerance.linear()));
    let hit_nodes = crate::merge::graph_split_node::split_edges_at_face_hits(
        &mut graph,
        &intersection_result.edge_face_hits,
        snap_tol,
        breps,
    );

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
    let rotated: Vec<crate::intersect::IntersectionRecord> = intersection_result
        .records
        .iter()
        .map(|r| rotate_closed_record_to_hits(r, &intersection_result.edge_face_hits, tolerance))
        .collect();
    let mut rotated = rotated;
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
    insert_intersection_curves(&mut graph, &rotated, snap_tol);
    let _ = &hit_nodes;
    // Split SHARED edges at record×record T-junctions: one record's
    // endpoint frequently lands mid-span of another (a dovetail
    // rail's flank trace ends ON the joint-shoulder trace), and the
    // walk cannot route the junction without a node-split there
    // (measured: the rail's floor strip never became a fragment and
    // the carve went non-manifold). Every node is a candidate — the
    // splitter skips nodes not on a shared edge's interior.
    let all_nodes: Vec<usize> = (0..graph.nodes.len()).collect();
    split_shared_edges_at_nodes(&mut graph, &all_nodes, snap_tol);
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
    if !intersection_result.edge_face_hits.is_empty()
        || intersection_result
            .coincidence
            .face_overlaps
            .iter()
            .any(|o| !o.full_overlap)
    {
        split_edges_at_shared_nodes(&mut graph, snap_tol, breps);
        if intersection_result
            .coincidence
            .face_overlaps
            .iter()
            .any(|o| !o.full_overlap)
        {
            remove_collinear_duplicates(&mut graph);
        }
    }

    // Sort radial edges at all nodes.
    sort_all_radial(&mut graph);

    // Walk loops to find face fragments.
    let overlap_pairs: Vec<_> = intersection_result
        .coincidence
        .face_overlaps
        .iter()
        .map(|o| (o.face_a, o.face_b))
        .collect();
    let mut walk_fallbacks: Vec<(crate::merge::graph::FaceRef, String)> = Vec::new();
    let fragments = walk_all_fragments_ext(
        &graph,
        original_a,
        original_b,
        !intersection_result.coincidence.face_overlaps.is_empty(),
        true,
        &overlap_pairs,
        &mut walk_fallbacks,
        tolerance,
    );

    // Phase 6 (CELLULAR): Classify fragments by cell membership.
    let selections = super::classify_fragments::classify_fragments_cellular(
        &fragments,
        &graph,
        original_a,
        original_b,
        op,
        tolerance,
        &intersection_result.coincidence.face_overlaps,
    )?;

    // Order-independent assembly from selected fragments.
    let kept: Vec<usize> = (0..fragments.len())
        .filter(|&i| selections[i].keep)
        .collect();

    let result = super::assemble::cellular_assemble(
        &graph,
        &fragments,
        &selections,
        original_a,
        original_b,
        *tolerance,
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

    Ok(result)
}
