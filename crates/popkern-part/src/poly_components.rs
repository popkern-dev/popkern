//! Solids of a polyhedron: one per edge-connected set of faces.

use std::collections::HashMap;

use popkern_core::topology::*;

/// One solid per edge-connected set of faces: the first keeps
/// `solid`/`shell`, each further set gets its own (a design body can be
/// several pieces; one shell over them breaks Euler-Poincaré).
pub(super) fn split_components(brep: &mut BRep, solid: SolidHandle, shell: ShellHandle) {
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
                    let twin = brep.half_edges[he].twin;
                    let other = brep.loops[brep.half_edges[twin].loop_].face;
                    if let Some(&j) = index.get(&other)
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
    let mut shells = vec![(solid, shell)];
    for _ in 1..count {
        let s = brep.solids.insert(Solid {
            shell: Default::default(),
        });
        let sh = brep.shells.insert(Shell {
            face: Default::default(),
            solid: s,
        });
        brep.solids[s].shell = sh;
        shells.push((s, sh));
    }
    let mut represented = vec![false; count];
    for (i, &f) in faces.iter().enumerate() {
        let c = component[i];
        let (_, sh) = shells[c];
        brep.faces[f].shell = sh;
        if !represented[c] {
            represented[c] = true;
            brep.shells[sh].face = f;
        }
    }
}
