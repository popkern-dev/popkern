//! Retirement of unglued construction-residual scaffolding.
//!
//! Apex construction (notably the V-free repeated-apex partition)
//! retires planning fragments — apex cells, displaced boundary
//! fragments, two-sided holding faces — into
//! [`TopologyRole::ConstructionResidual`] faces. Fragments that seal a
//! slit in the production shell (their sides twin production-owned
//! sides) are load-bearing membranes and must survive. Fragments with
//! no production-twinned side are pure scaffolding: they share at most
//! vertices with the production shell, and leaving them in the result
//! breaks downstream serving (STEP export refuses to synthesize a
//! plane for a two-vertex loop; tessellation and measure walk them as
//! geometry).
//!
//! This pass removes exactly the unglued fragments, transactionally:
//! the whole retirement rolls back if the cleaned model fails
//! validation. Gluing can only be broken by editing production
//! topology — which this pass never touches — so a single sweep is
//! complete.

use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle, TopologyRole};

/// Structured counts for the scaffolding retirement performed at the
/// end of a blend operation.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ResidualScaffoldingRetirementReport {
    /// Construction-residual faces removed (no production-twinned side).
    pub retired_faces: usize,
    /// Edges removed because both sides belonged to retired faces.
    pub retired_edges: usize,
    /// Vertices removed because no live half-edge referenced them.
    pub retired_vertices: usize,
    /// Construction-residual faces retained because at least one of
    /// their sides twins a production-owned side (slit membranes).
    pub retained_glued_faces: usize,
    /// True when the retirement was rolled back because the cleaned
    /// model failed validation; counts are zero in that case.
    pub rolled_back: bool,
}

/// Remove every [`TopologyRole::ConstructionResidual`] face none of
/// whose sides twin a production-owned side, together with the edges
/// and vertices that become unreferenced.
///
/// Invariants preserved: production faces, loops, half-edges, edges,
/// and their twin pairings are never modified; retained residual
/// membranes keep their loops intact. Sides of retained faces whose
/// twin is removed become boundary sides (their edge joins
/// `boundary_edges`), the same representation already used for
/// unowned residual sides. Postcondition: `BRep::validate()` passes on
/// the cleaned model, enforced by rollback.
pub(crate) fn retire_unglued_construction_residual_scaffolding(
    brep: &mut BRep,
) -> ResidualScaffoldingRetirementReport {
    let mut report = ResidualScaffoldingRetirementReport::default();
    // Side ownership is determined by the recorded `loop_` membership
    // of every half-edge, never by walking next-chains: a stale
    // residual face (one whose sides were re-owned by production
    // loops during a rewrite) has a chain that wanders live topology,
    // but it OWNS nothing — it is a ghost record and retires as a
    // face/loop pair without touching a single half-edge.
    let mut sides_by_loop: std::collections::BTreeMap<LoopHandle, Vec<HalfEdgeHandle>> =
        std::collections::BTreeMap::new();
    for (half_edge, half_edge_data) in &brep.half_edges {
        sides_by_loop
            .entry(half_edge_data.loop_)
            .or_default()
            .push(half_edge);
    }
    let mut retire_faces: Vec<(FaceHandle, Vec<LoopHandle>, Vec<HalfEdgeHandle>)> = Vec::new();
    for (face, face_data) in &brep.faces {
        if face_data.role != TopologyRole::ConstructionResidual {
            continue;
        }
        let mut loops = vec![face_data.outer_loop];
        loops.extend(face_data.inner_loops.iter().copied());
        let sides: Vec<HalfEdgeHandle> = loops
            .iter()
            .filter_map(|loop_handle| sides_by_loop.get(loop_handle))
            .flatten()
            .copied()
            .collect();
        let glued = sides.iter().any(|side| {
            let twin = brep.half_edges.get(*side).map(|data| data.twin);
            twin.and_then(|twin| brep.half_edges.get(twin))
                .and_then(|twin_data| brep.loops.get(twin_data.loop_))
                .and_then(|loop_data| brep.faces.get(loop_data.face))
                .is_some_and(|face_data| face_data.role != TopologyRole::ConstructionResidual)
        });
        if glued {
            report.retained_glued_faces += 1;
        } else {
            retire_faces.push((face, loops, sides));
        }
    }
    if retire_faces.is_empty() {
        return report;
    }

    let rollback = brep.clone();
    for (face, loops, sides) in &retire_faces {
        for side in sides {
            let Some(side_data) = brep.remove_half_edge_reanchored(*side) else {
                continue;
            };
            match brep.half_edges.get(side_data.twin) {
                Some(_) => {
                    // The twin survives (a retained membrane side or a
                    // side of a residual face retired later in this
                    // sweep); if it survives the sweep its edge is a
                    // boundary edge. The edge's canonical half-edge
                    // must not dangle on the removed side.
                    brep.boundary_edges.insert(side_data.edge);
                    if let Some(edge_data) = brep.edges.get_mut(side_data.edge)
                        && edge_data.half_edge == *side
                    {
                        edge_data.half_edge = side_data.twin;
                    }
                }
                None => {
                    if brep.edges.remove(side_data.edge).is_some() {
                        brep.boundary_edges.remove(&side_data.edge);
                        report.retired_edges += 1;
                    }
                }
            }
        }
        for loop_handle in loops {
            brep.loops.remove(*loop_handle);
        }
        brep.faces.remove(*face);
        report.retired_faces += 1;
    }

    // Drop boundary markers for edges that were fully removed after
    // their marker was set (twin removed later in the sweep).
    brep.boundary_edges
        .retain(|edge| brep.edges.contains_key(*edge));

    // Repair vertex outgoing pointers and drop unreferenced vertices.
    let mut outgoing_by_vertex = std::collections::BTreeMap::new();
    for (half_edge, half_edge_data) in &brep.half_edges {
        outgoing_by_vertex
            .entry(half_edge_data.origin)
            .or_insert(half_edge);
    }
    let vertices: Vec<_> = brep.vertices.keys().collect();
    for vertex in vertices {
        match outgoing_by_vertex.get(&vertex) {
            Some(outgoing) => {
                if !brep.half_edges.contains_key(brep.vertices[vertex].outgoing) {
                    brep.vertices[vertex].outgoing = *outgoing;
                }
            }
            None => {
                brep.vertices.remove(vertex);
                report.retired_vertices += 1;
            }
        }
    }

    // Repoint shells whose representative face was retired.
    let shells: Vec<_> = brep.shells.keys().collect();
    for shell in shells {
        if brep.faces.contains_key(brep.shells[shell].face) {
            continue;
        }
        if let Some(replacement) = brep
            .faces
            .iter()
            .find_map(|(face, data)| (data.shell == shell).then_some(face))
        {
            brep.shells[shell].face = replacement;
        }
    }

    if brep.validate().is_err() {
        *brep = rollback;
        return ResidualScaffoldingRetirementReport {
            rolled_back: true,
            ..Default::default()
        };
    }
    report
}
