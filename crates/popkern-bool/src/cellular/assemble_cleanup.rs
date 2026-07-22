//! Assembly cleanup and provenance: unwired-half-edge removal and
//! output-face provenance construction.
//!
//! Split from `cellular/assemble.rs` to keep it within the
//! structural line ratchet.

use std::collections::HashMap;

use popkern_core::constructors::builder::BRepBuilder;
use popkern_core::topology::{BRep, HalfEdgeHandle, VertexHandle};

use crate::merge::fragment_classify::FragmentSelection;
use crate::merge::graph::Operand;
use crate::merge::walk::Fragment;

/// Remove half-edges not wired into any loop, and their orphaned edges/vertices.
pub(super) fn cleanup_unwired(builder: &mut BRepBuilder) {
    let brep = builder.brep_mut();
    let default_loop = popkern_core::topology::LoopHandle::default();

    let dangling: Vec<HalfEdgeHandle> = brep
        .half_edges
        .iter()
        .filter(|(_, he)| he.loop_ == default_loop)
        .map(|(h, _)| h)
        .collect();

    let twins: Vec<(HalfEdgeHandle, HalfEdgeHandle)> = dangling
        .iter()
        .filter_map(|&h| brep.half_edges.get(h).map(|he| (h, he.twin)))
        .collect();

    for (h, twin_h) in &twins {
        if let Some(twin) = brep.half_edges.get_mut(*twin_h) {
            twin.twin = *twin_h;
        }
        // The surviving edge KEEPS its curve: edge geometry is a
        // property of the edge, not of its twin wiring — clearing it
        // here silently degraded seam-style survivors (a boss/plate
        // intersection circle lost its geometry entirely).
        brep.half_edges.remove(*h);
    }

    // Remove fully-orphaned edges
    let dead_edges: Vec<_> = brep
        .edges
        .iter()
        .filter(|(_, e)| !brep.half_edges.contains_key(e.half_edge))
        .filter(|(eh, _)| !brep.half_edges.iter().any(|(_, he)| he.edge == *eh))
        .map(|(eh, _)| eh)
        .collect();
    for eh in dead_edges {
        brep.edges.remove(eh);
    }

    // Fix vertex outgoing handles
    let outgoing: HashMap<VertexHandle, HalfEdgeHandle> = brep
        .half_edges
        .iter()
        .map(|(h, he)| (he.origin, h))
        .collect();
    for (vh, v) in &mut brep.vertices {
        if !brep.half_edges.contains_key(v.outgoing) {
            v.outgoing = outgoing
                .get(&vh)
                .copied()
                .unwrap_or(HalfEdgeHandle::default());
        }
    }

    // Remove orphan vertices
    let orphans: Vec<VertexHandle> = brep
        .vertices
        .iter()
        .filter(|(_, v)| v.outgoing == HalfEdgeHandle::default())
        .map(|(h, _)| h)
        .collect();
    for vh in orphans {
        brep.vertices.remove(vh);
    }
}

/// Build entity provenance map from assembly information.
///
/// Records which input face each output face came from, and which
/// graph edge source each output vertex came from.
/// Build entity provenance from assembly data.
///
/// Populates face provenance by matching output faces to their source
/// fragments. Called after assembly when provenance tracking is needed.
pub(crate) fn build_face_provenance(
    result: &BRep,
    kept: &[usize],
    fragments: &[Fragment],
    _selections: &[FragmentSelection],
) -> crate::provenance::EntityMap {
    use crate::provenance::types::{EntityMap, FaceProvenance, SourceBody};

    let mut map = EntityMap::default();

    let result_faces: Vec<_> = result.faces.keys().collect();

    // Match output faces to input faces via fragment ordering.
    // The assembly creates faces in the same order as `kept` iteration.
    for (i, &fi) in kept.iter().enumerate() {
        if i >= result_faces.len() {
            break;
        }
        let frag = &fragments[fi];
        let source = match frag.face.operand {
            Operand::A => SourceBody::A,
            Operand::B => SourceBody::B,
        };
        map.faces.push((
            result_faces[i],
            FaceProvenance::Kept {
                source,
                original_face: frag.face.face,
            },
        ));
    }

    map
}
