//! Assembly cleanup and provenance: unwired-half-edge removal and
//! output-face provenance construction.
//!
//! Split from `cellular/assemble.rs` to keep it within the
//! structural line ratchet.

use popkern_math::hash::FastHashMap as HashMap;

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
        .filter_map(|&h| {
            brep.half_edges
                .get(h)
                .and_then(|he| he.mate())
                .map(|t| (h, t))
        })
        .collect();

    for (h, twin_h) in &twins {
        if let Some(twin) = brep.half_edges.get_mut(*twin_h) {
            twin.set_mate(*twin_h);
        }
        // The surviving edge KEEPS its curve: edge geometry is a
        // property of the edge, not of its twin wiring — clearing it
        // here silently degraded seam-style survivors (a boss/plate
        // intersection circle lost its geometry entirely).
        brep.remove_half_edge_reanchored(*h);
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
/// Every component becomes its own solid. A cavity (a component whose
/// faces enclose negative volume) is nested in the solid holding it after
/// the output repairs, on the finished body (`postprocess::cavity_nesting`,
/// the one nesting authority, PK-1293); the mid-assembly ray parity this
/// replaced read unfinished trims (PK-217's rule, retired).
pub(super) fn partition_shells(builder: &mut BRepBuilder) -> bool {
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
        let Some(twin) = he
            .mate()
            .filter(|t| *t != h && brep.half_edges.contains_key(*t))
        else {
            continue;
        };
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
    let mut component_of: HashMap<usize, usize> = HashMap::default();
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

    let face_component: Vec<(popkern_core::topology::FaceHandle, usize)> = faces
        .iter()
        .enumerate()
        .map(|(i, &f)| (f, component_of[&find(&mut parent, i)]))
        .collect();

    // Every component is its own solid: nesting has one authority, the
    // post-repair cavity nesting (`postprocess::cavity_nesting`), which
    // reads the finished body. Mid-assembly parity on unfinished trims
    // declined whenever a nested component offered no interior point
    // and served two solids as one shell (a ball poking out through a
    // cylinder's cap, minus the cylinder: corpus scene 1337), and read a
    // void cut from two balls by a third Outside them (PK-1293).
    let encloser: Vec<Option<usize>> = vec![None; roots.len()];
    // Wire shells: every component gets a shell; solid-owning
    // components (even enclosure depth) get their own solid, cavity
    // components join their encloser's solid. The first solid-owning
    // component reuses the existing shell/solid.
    let first_own = (0..roots.len()).find(|&c| encloser[c].is_none());
    let Some(first_own) = first_own else {
        return false; // every component claims to be a cavity — lie
    };
    let Some(existing_shell) = brep.face_back_shell(faces[0]) else {
        return false;
    };
    let Some(existing_solid) = brep.shell_region(existing_shell) else {
        return false;
    };
    let Some(existing_sides) = brep.face_sides(faces[0]) else {
        return false;
    };
    let mut shells: Vec<Option<popkern_core::topology::FaceSides>> = vec![None; roots.len()];
    let mut solids: Vec<Option<popkern_core::topology::RegionHandle>> = vec![None; roots.len()];
    shells[first_own] = Some(existing_sides);
    solids[first_own] = Some(existing_solid);
    for c in 0..roots.len() {
        if c == first_own || encloser[c].is_some() {
            continue;
        }
        let solid = builder.add_solid();
        solids[c] = Some(solid);
        let shell = builder.add_shell(solid);
        shells[c] = Some(builder.shell_sides(shell));
    }
    let mut cavity_count = 0usize;
    for c in 0..roots.len() {
        if let Some(e) = encloser[c] {
            let Some(owner_solid) = solids[e] else {
                return false;
            };
            let shell = builder.add_shell(owner_solid);
            shells[c] = Some(builder.shell_sides(shell));
            cavity_count += 1;
        }
    }
    // A region's first shell stays its outer shell: cavity shells
    // added after it do not usurp the slot.
    let brep = builder.brep_mut();
    // Re-anchor every shell's representative face (component 0's may
    // have pointed at a face now moved to another shell).
    let mut anchored = vec![false; roots.len()];
    for (f, c) in face_component {
        let Some(sides) = shells[c] else { return false };
        brep.set_face_sides(f, sides);
        if !anchored[c] {
            brep.set_sides_anchor(sides, f);
            anchored[c] = true;
        }
    }
    if popkern_core::trace::trace_enabled("boolean") {
        eprintln!(
            "  [assemble] shell partition: {} component(s) → {} solid(s), {} cavity shell(s)",
            roots.len(),
            roots.len() - cavity_count,
            cavity_count
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
            _ /* B */ => SourceBody::B,
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
