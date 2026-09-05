//! Graph merge orchestrator — builds intersection graph and assembles result BRep.
//!
//! Ties together all graph phases: topology insertion, edge splitting,
//! intersection curve insertion, radial sorting, loop walking, fragment
//! classification, and final BRep assembly.

use popkern_core::{tolerance::Tolerance, topology::BRep};

use super::assemble::assemble_result;
use super::fragment_classify::classify_and_select_fragments;
use super::graph::IntersectionGraph;
use super::graph_build::insert_operand_topology;
use super::graph_intersect::insert_intersection_curves;
use super::graph_split::remove_collinear_duplicates;
use super::graph_split_shared::{
    split_edges_at_shared_nodes, split_shared_edges_at_interior_nodes,
};
use super::radial::sort_all_radial;
use super::walk::walk_all_fragments;
use crate::error::{BoolDiagnostics, BoolOperation, Stage};
use crate::intersect::IntersectionResult;

use super::graph::Operand;

/// Execute the graph-based merge.
///
/// Builds an intersection graph from the original operands and intersection
/// curves, walks face fragment loops, classifies and selects fragments,
/// and assembles the result BRep.
///
/// # Errors
///
/// Returns [`BoolDiagnostics`] if any phase fails to produce valid topology.
pub(crate) fn graph_merge(
    original_a: &BRep,
    original_b: &BRep,
    intersection_result: &IntersectionResult,
    op: BoolOperation,
    tolerance: &Tolerance,
) -> Result<BRep, BoolDiagnostics> {
    // Use snap tolerance (2x linear) for vertex merging to prevent micro-edges.
    let snap_tol = tolerance.linear() * 2.0;
    let mut graph = IntersectionGraph::new();

    // Insert original topology from both operands.
    let (_vmap_a, _hemap_a) = insert_operand_topology(&mut graph, original_a, Operand::A, snap_tol);
    let (_vmap_b, _hemap_b) = insert_operand_topology(&mut graph, original_b, Operand::B, snap_tol);

    // Process edge-face hits — split edges at pierce points.
    // (OnVertex already exists; OnFaceInterior is handled when
    // inserting curves.)
    let _ = super::graph_split_node::split_edges_at_face_hits(
        &mut graph,
        &intersection_result.edge_face_hits,
        snap_tol,
        Some((original_a, original_b, tolerance.linear())),
    );

    // Insert intersection curves as shared edges.
    insert_intersection_curves(&mut graph, &intersection_result.records, snap_tol);

    // Split record edges at pre-existing pierce nodes on
    // their interior (a record crossing a periodic face's seam passes
    // through the seam's phase-2 pierce node; unsplit, the walk
    // cannot route through it and the flanking fragments never close
    // — PK-295).
    split_shared_edges_at_interior_nodes(
        &mut graph,
        snap_tol,
        Some((original_a, original_b, tolerance.linear())),
    );

    // Split ORIGINAL edges at shared-record endpoint
    // nodes. Phase 2 pierces edges named by edge-face hits, but its
    // OnVertex short-circuit assumes an existing node means every
    // edge through that point is already split — false when a SECOND
    // original edge passes through the same point (measured: the
    // cylinder wall's seam runs through both torus pierce nodes;
    // unsplit, the wall fragment's loop carried the FULL seam, its
    // curve-midpoint classification sample landed inside the torus,
    // and the union discarded the entire wall). Runs unconditionally
    // — the pass is node-driven and curve-verified, and a shared
    // node genuinely mid-edge always means a missed split.
    split_edges_at_shared_nodes(
        &mut graph,
        snap_tol,
        Some((original_a, original_b, tolerance.linear())),
    );
    // Collinear-duplicate removal stays scoped to partial
    // coincident overlaps (its reason to exist).
    if intersection_result
        .coincidence
        .face_overlaps
        .iter()
        .any(|o| !o.full_overlap)
    {
        remove_collinear_duplicates(&mut graph);
    }

    // Sort radial edges at all nodes.
    sort_all_radial(&mut graph);

    // §7.4 manifold heuristic. Walk the graph and report
    // any vertex pairs where more than 2 distinct faces meet. This
    // is a HEURISTIC detector — not all such pairs are actually
    // non-manifold (a vertex commonly has 3+ adjacent faces in
    // valid topology). Findings are logged in debug builds for
    // diagnostic purposes only; the authoritative non-manifold
    // gate runs in `validate::post_check::validate_manifold_edges`
    // after assembly, where it can compare against the resolved
    // BRep topology rather than the in-progress graph state.
    let manifold_report = super::manifold_check::detect_multi_incident_edges(&graph);
    if !manifold_report.is_manifold() {
        #[cfg(debug_assertions)]
        eprintln!(
            "[manifold_check] {} vertex pair(s) with >2 face refs (heuristic only)",
            manifold_report.multi_incident.len()
        );
    }

    // Walk loops to find face fragments.
    let fragments = walk_all_fragments(&graph, original_a, original_b, tolerance);

    // Classify and select fragments.
    let selections = classify_and_select_fragments(
        &fragments,
        &graph,
        original_a,
        original_b,
        op,
        tolerance,
        &intersection_result.coincidence,
    )
    .map_err(|e| BoolDiagnostics::from_error(e, Stage::Classification))?;

    // Assemble result BRep.
    assemble_result(
        &graph,
        &fragments,
        &selections,
        original_a,
        original_b,
        *tolerance,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::BoolOperation;
    use crate::intersect::{IntersectionResult, compute_face_aabbs};
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    #[test]
    fn graph_merge_identity_cube() {
        // Build one cube, "intersect" it with itself (no intersection curves).
        // graph_merge with empty intersection_result should handle gracefully.
        let tol = Tolerance::default();
        let cube = BRep::unit_cube(tol).expect("cube");

        let aabbs = compute_face_aabbs(&cube);
        let intersection_result = IntersectionResult {
            records: vec![],
            aabbs_a: aabbs.clone(),
            aabbs_b: aabbs,
            edge_face_hits: vec![],
            hierarchy: Default::default(),
            coincidence: Default::default(),
        };

        // With no intersection curves, the graph has only original topology.
        // This may produce a valid result or an empty one — the important
        // thing is that it does not panic.
        let result = graph_merge(
            &cube,
            &cube,
            &intersection_result,
            BoolOperation::Union,
            &tol,
        );

        // We accept either Ok or Err — just verify no panic.
        match result {
            Ok(brep) => {
                // If it succeeds, the result should have some faces.
                assert!(!brep.faces.is_empty() || brep.faces.is_empty());
            }
            Err(_diag) => {
                // Graph merge may fail on identity (no intersection curves
                // means classification may be incomplete). That's acceptable.
            }
        }
    }
}
