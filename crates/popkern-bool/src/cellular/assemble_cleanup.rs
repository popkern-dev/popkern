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

/// Partition the assembled faces into edge-connected components and
/// register one shell (and solid) per component.
///
/// A boolean can DISCONNECT the workpiece — a flush through-slot
/// subtract leaves two closed boxes — and an assembly that files
/// every face under one shell fails Euler-Poincaré by construction
/// (χ = 2 per closed component, the formula expects 2 per SHELL), so
/// the whole cellular result fell back to graph-merge (measured on
/// the slotted-block corpus part: correct 12-face fragment population
/// discarded, graph-merge served 48 open edges).
///
/// Components here are mutually EXTERIOR pieces: each becomes its own
/// solid. A component whose AABB nests strictly inside another's is
/// ambiguous between a cavity shell (same solid) and a floating piece
/// in a through-hole (own solid) — deciding needs a containment
/// oracle, and misfiling changes measured mass properties. That
/// configuration DECLINES the partition (returns false, single-shell
/// result stands, the Euler gate routes it to graph-merge exactly as
/// before) with the reason observable under PK_BOOL_TRACE.
pub(super) fn partition_shells(builder: &mut BRepBuilder, tolerance_linear: f64) -> bool {
    let brep = builder.brep_mut();
    let faces: Vec<_> = brep.faces.keys().collect();
    if faces.len() < 2 {
        return true;
    }
    let index: HashMap<_, usize> = faces.iter().enumerate().map(|(i, &f)| (f, i)).collect();

    // Union-find over faces joined by a shared (twinned) edge.
    let mut parent: Vec<usize> = (0..faces.len()).collect();
    fn find(parent: &mut [usize], i: usize) -> usize {
        let mut r = i;
        while parent[r] != r {
            r = parent[r];
        }
        let mut c = i;
        while parent[c] != r {
            let next = parent[c];
            parent[c] = r;
            c = next;
        }
        r
    }
    for (h, he) in brep.half_edges.iter() {
        let twin = he.twin;
        if twin == h || !brep.half_edges.contains_key(twin) {
            continue;
        }
        let (Some(l1), Some(l2)) = (
            brep.loops.get(he.loop_),
            brep.loops.get(brep.half_edges[twin].loop_),
        ) else {
            continue;
        };
        if let (Some(&i1), Some(&i2)) = (index.get(&l1.face), index.get(&l2.face)) {
            let (r1, r2) = (find(&mut parent, i1), find(&mut parent, i2));
            if r1 != r2 {
                parent[r1] = r2;
            }
        }
    }
    let mut component_of: HashMap<usize, usize> = HashMap::new();
    let mut roots: Vec<usize> = Vec::new();
    for i in 0..faces.len() {
        let r = find(&mut parent, i);
        let c = *component_of.entry(r).or_insert_with(|| {
            roots.push(r);
            roots.len() - 1
        });
        let _ = c;
    }
    if roots.len() < 2 {
        return true;
    }

    // Component AABBs from wired half-edge origins.
    let mut aabbs: Vec<([f64; 3], [f64; 3])> = vec![([f64::MAX; 3], [f64::MIN; 3]); roots.len()];
    for (_, he) in brep.half_edges.iter() {
        let Some(lp) = brep.loops.get(he.loop_) else {
            continue;
        };
        let Some(&fi) = index.get(&lp.face) else {
            continue;
        };
        let r = find(&mut parent, fi);
        let c = component_of[&r];
        let p = brep.vertices[he.origin].position;
        for (k, v) in [p.x, p.y, p.z].into_iter().enumerate() {
            aabbs[c].0[k] = aabbs[c].0[k].min(v);
            aabbs[c].1[k] = aabbs[c].1[k].max(v);
        }
    }
    let tol = tolerance_linear.max(1e-12);
    let nested = |a: &([f64; 3], [f64; 3]), b: &([f64; 3], [f64; 3])| -> bool {
        (0..3).all(|k| a.0[k] > b.0[k] - tol && a.1[k] < b.1[k] + tol)
    };
    for i in 0..roots.len() {
        for j in 0..roots.len() {
            if i != j && nested(&aabbs[i], &aabbs[j]) {
                if std::env::var_os("PK_BOOL_TRACE").is_some() {
                    eprintln!(
                        "  [assemble] shell partition DECLINED: component {i} AABB nests in {j} (cavity vs floating piece needs containment)"
                    );
                }
                return false;
            }
        }
    }

    // Reassign: component 0 keeps the existing shell; each further
    // component gets a fresh solid + shell.
    let face_component: Vec<(popkern_core::topology::FaceHandle, usize)> = faces
        .iter()
        .enumerate()
        .map(|(i, &f)| (f, component_of[&find(&mut parent, i)]))
        .collect();
    let existing_shell = brep.faces[faces[0]].shell;
    let mut shells: Vec<popkern_core::topology::ShellHandle> = vec![existing_shell];
    for _ in 1..roots.len() {
        let solid = builder.add_solid();
        shells.push(builder.add_shell(solid));
    }
    let brep = builder.brep_mut();
    // Re-anchor every shell's representative face (component 0's may
    // have pointed at a face now moved to another shell).
    let mut anchored = vec![false; shells.len()];
    for (f, c) in face_component {
        brep.faces[f].shell = shells[c];
        if !anchored[c] {
            brep.shells[shells[c]].face = f;
            anchored[c] = true;
        }
    }
    if std::env::var_os("PK_BOOL_TRACE").is_some() {
        eprintln!(
            "  [assemble] shell partition: {} disconnected component(s) → one shell/solid each",
            shells.len()
        );
    }
    true
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
