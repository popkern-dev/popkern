//! Solids of a polyhedron: one per edge-connected set of faces.

use std::collections::HashMap;

use popkern_core::topology::*;

/// One solid per edge-connected set of faces: the first keeps
/// `sides`, each further set gets its own (a design body can be
/// several pieces; one shell over them breaks Euler-Poincaré).
pub(super) fn split_components(brep: &mut BRep, first: FaceSides) {
    let faces: Vec<FaceHandle> = brep.faces.keys().collect();
    let index: HashMap<FaceHandle, usize> =
        faces.iter().enumerate().map(|(i, &f)| (f, i)).collect();
    let mut component = vec![usize::MAX; faces.len()];
    let mut count = 0;
    for start in 0..faces.len() {
        if component[start] != usize::MAX {
            continue;
        }
        component[start] = count;
        let mut stack = vec![start];
        while let Some(i) = stack.pop() {
            let face = &brep.faces[faces[i]];
            let loops: Vec<LoopHandle> = std::iter::once(face.outer_loop)
                .chain(face.inner_loops.iter().copied())
                .collect();
            for lh in loops {
                let first = brep.loops[lh].half_edge;
                let mut he = first;
                loop {
                    let other = brep
                        .mate(he)
                        .map(|twin| brep.loops[brep.half_edges[twin].loop_].face);
                    if let Some(&j) = other.and_then(|other| index.get(&other))
                        && component[j] == usize::MAX
                    {
                        component[j] = count;
                        stack.push(j);
                    }
                    he = brep.half_edges[he].next;
                    if he == first {
                        break;
                    }
                }
            }
        }
        count += 1;
    }
    let mut sides = vec![first];
    for _ in 1..count {
        sides.push(brep.add_solid_region(Default::default()).1);
    }
    let mut represented = vec![false; count];
    for (i, &f) in faces.iter().enumerate() {
        let c = component[i];
        brep.set_face_sides(f, sides[c]);
        if !represented[c] {
            represented[c] = true;
            brep.set_sides_anchor(sides[c], f);
        }
    }
}
