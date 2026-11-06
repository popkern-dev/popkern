//! Building one cell complex over N bodies.
//!
//! The two-operand pipeline's stages, run once over all bodies: every
//! body's topology under its own operand index; for each body pair
//! whose boxes overlap, its intersection records (closed records closed
//! at an edge-face hit on them) with the pair's operands, and its
//! operand edges split at the pair's hits; nodes where records of
//! different pairs cross on a common face; the graph's own junction
//! splits; one chart walk; one complex through membership. Vertex
//! identity is the graph's snap-then-adopt identity in this stage (the
//! vertex table resolves one body pair).

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::{BodySet, RegionInfo};
use crate::complex::model::{CellComplex, SideIndex};
use crate::error::BoolError;
use crate::merge::graph::{FacePair, FaceRef, IntersectionGraph, Operand};

/// Whether two face-box lists overlap anywhere.
fn boxes_overlap(
    a: &[(
        popkern_core::topology::FaceHandle,
        popkern_core::geometry::intersection::surface_subdivision::Aabb,
    )],
    b: &[(
        popkern_core::topology::FaceHandle,
        popkern_core::geometry::intersection::surface_subdivision::Aabb,
    )],
    reach: f64,
) -> bool {
    a.iter()
        .any(|(_, x)| b.iter().any(|(_, y)| x.overlaps(y, reach)))
}

fn declined(stage: &str, reason: impl std::fmt::Display) -> BoolError {
    BoolError::MergeFailed {
        stage: stage.to_string(),
        reason: format!("partition: {reason}"),
    }
}

/// What a partition keeps to answer requests: the graph, its
/// fragments, the complex through membership, and the regions.
pub(super) struct Built {
    pub graph: IntersectionGraph,
    pub fragments: Vec<crate::merge::walk::Fragment>,
    pub complex: crate::complex::pipeline::ComplexRegions,
    pub regions: Vec<RegionInfo>,
    /// The complex shells bounding each region, in `regions` order,
    /// the outer shell first (`region_nesting`).
    pub region_shells: Vec<Vec<usize>>,
}

/// Build the complex of `bodies` and read off its bounded regions, each
/// with its membership and an interior point.
/// `lamina_edges[k]` are operand `k`'s open boundary edges when it is a
/// sheet entering as a lamina (empty otherwise).
pub(super) fn build(
    bodies: &[&BRep],
    lamina_edges: &[std::collections::BTreeSet<popkern_core::topology::EdgeHandle>],
    tolerance: &Tolerance,
) -> Result<Built, BoolError> {
    let linear = tolerance.linear();
    let snap = linear * 2.0;
    let boxes: Vec<_> = bodies
        .iter()
        .map(|b| crate::intersect::compute_face_aabbs(b))
        .collect();
    let mut graph = IntersectionGraph::new();
    for (i, body) in bodies.iter().enumerate() {
        crate::merge::graph_build::insert_operand_topology_scoped(
            &mut graph,
            body,
            Operand::from_index(i),
            snap,
            None,
            None,
        );
    }
    let breps = Some((bodies, linear));
    // Coincident faces of every pair, and whether any only partly
    // overlap (their records then run along the faces' own edges).
    let mut overlaps: Vec<FacePair> = Vec::new();
    let mut partial = false;
    for i in 0..bodies.len() {
        for j in i + 1..bodies.len() {
            if !boxes_overlap(&boxes[i], &boxes[j], linear) {
                continue;
            }
            let ix = crate::intersect::find_intersections(bodies[i], bodies[j], tolerance)
                .map_err(|d| d.into_error())?;
            let (oi, oj) = (Operand::from_index(i), Operand::from_index(j));
            let face = |operand, face| FaceRef { operand, face };
            for o in &ix.coincidence.face_overlaps {
                overlaps.push((face(oi, o.face_a), face(oj, o.face_b)));
                partial |= !o.full_overlap;
            }
            let mut records: Vec<_> = ix
                .records
                .iter()
                .map(|r| {
                    crate::cellular::record_rotation::rotate_closed_record_to_hits(
                        r,
                        &ix.edge_face_hits,
                        tolerance,
                    )
                })
                .collect();
            // A pair with coincident faces drops its trim residues, as
            // the two-operand pipeline does: the coincidence's own
            // records bound the shared patch.
            if !ix.coincidence.face_overlaps.is_empty() {
                records.retain(|r| r.residue_owner_a.is_none());
            }
            crate::merge::graph_split_node::split_edges_at_face_hits(
                &mut graph,
                &ix.edge_face_hits,
                snap,
                breps,
            );
            crate::merge::graph_intersect::insert_pair_curves(
                &mut graph,
                &records,
                snap,
                None,
                [Operand::from_index(i), Operand::from_index(j)],
            );
        }
    }
    super::crossings::split_at_record_crossings(&mut graph, snap, tolerance);
    let all_nodes: Vec<usize> = (0..graph.nodes.len()).collect();
    crate::merge::graph_split::split_shared_edges_at_nodes(&mut graph, &all_nodes, snap);
    crate::merge::graph_split_shared::split_edges_at_shared_nodes(&mut graph, snap, breps);
    if partial {
        crate::merge::graph_duplicates::remove_collinear_duplicates(&mut graph, linear);
    }
    crate::merge::radial::sort_all_radial(&mut graph);
    mark_lamina_boundary_nodes(&mut graph, lamina_edges);

    let mut fallbacks = Vec::new();
    let fragments = crate::merge::walk::walk_all_fragments_ext(
        &graph,
        bodies,
        !overlaps.is_empty(),
        true,
        &overlaps,
        &mut fallbacks,
        None,
        tolerance,
    );
    let fragments =
        crate::complex::hole_phantoms::rejoin_own_holes(&graph, &fragments, bodies, linear)
            .unwrap_or(fragments);
    let box_slices: Vec<&[_]> = boxes.iter().map(|b| &b[..]).collect();
    let frontier = std::collections::BTreeSet::new();
    let laminae: Vec<bool> = lamina_edges.iter().map(|e| !e.is_empty()).collect();
    let input = crate::complex::pipeline::ComplexInput {
        graph: &graph,
        fragments: &fragments,
        operands: bodies,
        boxes: &box_slices,
        overlaps: &overlaps,
        tolerance,
        frontier: &frontier,
        laminae: &laminae,
    };
    let complex = crate::complex::pipeline::complex_regions(&input, "partition")
        .map_err(|d| declined("Complex", format!("{d:?}")))?;
    let c = &complex.lift.complex;
    let mut out = Vec::new();
    let mut region_shells = Vec::new();
    for (si, shell) in c.shells.iter().enumerate() {
        let Some(membership) = shell.membership.as_ref() else {
            return Err(declined(
                "Complex",
                format!("region {si} has no membership"),
            ));
        };
        let bodies_in: Vec<usize> = membership.iter().collect();
        if bodies_in.is_empty() {
            continue;
        }
        let interior_point = interior_point(c, bodies, linear, si)
            .ok_or_else(|| declined("Complex", format!("region {si} has no interior point")))?;
        out.push(RegionInfo {
            bodies: BodySet::of(&bodies_in),
            interior_point,
        });
        region_shells.push(vec![si]);
    }
    Ok(Built {
        graph,
        fragments,
        complex,
        regions: out,
        region_shells,
    })
}

/// A point strictly inside region `shell`: a point just inside one of
/// its faces (as membership probes it), stepped off the face into the
/// region.
fn interior_point(
    c: &CellComplex,
    operands: &[&BRep],
    linear: f64,
    shell: usize,
) -> Option<popkern_core::geometry::Point3> {
    for &(face, side) in c.shells.get(shell)?.sides.iter().take(8) {
        let f = &c.faces[face.index()];
        let owner = f.owners.first()?;
        let surf = operands
            .get(owner.operand)?
            .faces
            .get(owner.face)?
            .surface
            .as_ref()?;
        let Some(at) = crate::complex::same_domain_check::inner_point(c, operands, linear, face)
        else {
            continue;
        };
        let n = surf.folded_normal(at.uv.0, at.uv.1)?;
        let outward = if f.folded_is_outward { n } else { n * -1.0 };
        let toward = if side == SideIndex::Ahead {
            outward
        } else {
            outward * -1.0
        };
        return Some(at.point + toward * (at.fin_length * 1e-3).max(linear * 50.0));
    }
    None
}

/// Record every node on a sheet operand's free boundary edges (PK-978).
fn mark_lamina_boundary_nodes(
    graph: &mut IntersectionGraph,
    lamina_edges: &[std::collections::BTreeSet<popkern_core::topology::EdgeHandle>],
) {
    use crate::merge::graph::EdgeSource;
    if lamina_edges.iter().all(|e| e.is_empty()) {
        return;
    }
    let mut marked = Vec::new();
    for he in &graph.half_edges {
        let (operand, edge) = match he.source {
            EdgeSource::Original(op, e)
            | EdgeSource::Split {
                operand: op,
                original: e,
            } => (op, e),
            EdgeSource::Shared { .. } => continue,
        };
        let k = operand.index();
        if lamina_edges
            .get(k)
            .is_some_and(|edges| edges.contains(&edge))
        {
            marked.push(he.origin);
            marked.push(graph.half_edges[he.twin].origin);
        }
    }
    graph.lamina_boundary_nodes.extend(marked);
}
