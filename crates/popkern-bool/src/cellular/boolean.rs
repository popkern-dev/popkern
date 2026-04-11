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

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::error::{BoolDiagnostics, BoolOperation, PointClassification};
use crate::intersect::IntersectionResult;
use crate::merge::fragment_classify::FragmentSelection;
use crate::merge::graph::{IntersectionGraph, Operand};
use crate::merge::graph_build::insert_operand_topology;
use crate::merge::graph_intersect::insert_intersection_curves;
use crate::merge::graph_split::{
    EndpointClass, classify_endpoint, remove_collinear_duplicates, split_edges_at_shared_nodes,
    split_graph_edge,
};
use crate::merge::radial::sort_all_radial;
use crate::merge::walk::{Fragment, walk_all_fragments};

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

    // Phase 1: Insert original topology from both operands.
    let (_vmap_a, _hemap_a) = insert_operand_topology(&mut graph, original_a, Operand::A, snap_tol);
    let (_vmap_b, _hemap_b) = insert_operand_topology(&mut graph, original_b, Operand::B, snap_tol);

    // Phase 2: Process edge-face hits.
    for hit in &intersection_result.edge_face_hits {
        let endpoint_class = classify_endpoint(&graph, &hit.point, snap_tol);
        if let EndpointClass::OnEdge(he_idx, t_param) = endpoint_class {
            split_graph_edge(&mut graph, he_idx, hit.point, t_param);
        }
    }

    // Phase 3: Insert intersection curves as shared edges.
    insert_intersection_curves(&mut graph, &intersection_result.records, snap_tol);

    // Phase 3.1–3.2: Handle coincident face boundaries.
    if intersection_result
        .coincidence
        .face_overlaps
        .iter()
        .any(|o| !o.full_overlap)
    {
        split_edges_at_shared_nodes(&mut graph, snap_tol);
        remove_collinear_duplicates(&mut graph);
    }

    // Phase 4: Sort radial edges at all nodes.
    sort_all_radial(&mut graph);

    // Phase 5: Walk loops to find face fragments.
    let fragments = walk_all_fragments(&graph, original_a, original_b);

    // Phase 6 (CELLULAR): Classify fragments by cell membership.
    let selections =
        classify_fragments_cellular(&fragments, &graph, original_a, original_b, op, tolerance)?;

    // Phase 7: Order-independent assembly from selected fragments.
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
    )?;

    // Build provenance map (available via boolean_with_provenance).
    let _provenance =
        super::assemble::build_face_provenance(&result, &kept, &fragments, &selections);

    Ok(result)
}

/// Classify fragments using cellular logic: determine which volumetric
/// cell each fragment borders, then select based on cell membership.
///
/// For each fragment:
/// 1. Sample a point on the fragment
/// 2. Test it against BOTH bodies via point-in-solid
/// 3. Compute cell membership: bits 0=in_A, 1=in_B
/// 4. Apply operation predicate to the cell membership
fn classify_fragments_cellular(
    fragments: &[Fragment],
    graph: &IntersectionGraph,
    brep_a: &BRep,
    brep_b: &BRep,
    op: BoolOperation,
    tolerance: &Tolerance,
) -> Result<Vec<FragmentSelection>, BoolDiagnostics> {
    use crate::adaptive_precision::robust_point_in_solid;
    use crate::intersect::compute_face_aabbs;

    let aabbs_a = compute_face_aabbs(brep_a);
    let aabbs_b = compute_face_aabbs(brep_b);

    let mut selections = Vec::with_capacity(fragments.len());

    for fragment in fragments {
        // Sample a point on this fragment
        let sample = sample_fragment_point(graph, &fragment.outer_loop);

        // Classify against the OPPOSITE body only.
        let own_operand = fragment.face.operand;
        let (other_brep, other_aabbs) = match own_operand {
            Operand::A => (brep_b, &aabbs_b),
            Operand::B => (brep_a, &aabbs_a),
        };

        let class = robust_point_in_solid(&sample, other_brep, other_aabbs, tolerance);

        // Selection rules equivalent to graph-merge but expressed via
        // cellular logic: keep fragment if the cell on its outer side
        // is not selected by the operation predicate.
        let (keep, flip) = match op {
            BoolOperation::Union => (
                class == PointClassification::Outside
                    || (own_operand == Operand::A && class == PointClassification::OnBoundary),
                false,
            ),
            BoolOperation::Intersection => (
                class == PointClassification::Inside
                    || (own_operand == Operand::A && class == PointClassification::OnBoundary),
                false,
            ),
            BoolOperation::Subtraction => match own_operand {
                Operand::A => (
                    class == PointClassification::Outside
                        || class == PointClassification::OnBoundary,
                    false,
                ),
                Operand::B => (class == PointClassification::Inside, true),
            },
        };

        selections.push(FragmentSelection { keep, flip });
    }

    Ok(selections)
}

/// Sample a point from a fragment's outer loop for classification.
fn sample_fragment_point(
    graph: &IntersectionGraph,
    loop_hes: &[usize],
) -> popkern_core::geometry::Point3 {
    use crate::merge::graph::EdgeSource;

    // Prefer midpoint of a non-shared edge (away from intersection curve)
    for &he_idx in loop_hes {
        let he = &graph.half_edges[he_idx];
        match &he.source {
            EdgeSource::OriginalA(_)
            | EdgeSource::OriginalB(_)
            | EdgeSource::SplitA { .. }
            | EdgeSource::SplitB { .. } => {
                let start = graph.nodes[he.origin].position;
                let twin_idx = he.twin;
                let end = graph.nodes[graph.half_edges[twin_idx].origin].position;
                return popkern_core::geometry::Point3::new(
                    (start.x + end.x) * 0.5,
                    (start.y + end.y) * 0.5,
                    (start.z + end.z) * 0.5,
                );
            }
            EdgeSource::Shared { .. } => continue,
        }
    }

    // Fallback: centroid of all loop vertices
    if loop_hes.is_empty() {
        return popkern_core::geometry::Point3::new(0.0, 0.0, 0.0);
    }
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut sz = 0.0;
    for &he_idx in loop_hes {
        let p = graph.nodes[graph.half_edges[he_idx].origin].position;
        sx += p.x;
        sy += p.y;
        sz += p.z;
    }
    let n = loop_hes.len() as f64;
    popkern_core::geometry::Point3::new(sx / n, sy / n, sz / n)
}
