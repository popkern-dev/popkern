//! FaceInfo (plan §4.3): what the resolution put on each operand face,
//! the input to imprint.
//!
//! Per face: the vertices **on** its boundary (a member on one of its
//! edges or vertices), the vertices **in** it (a member on the face and
//! none on its boundary), the pave blocks of its own edges the
//! resolution cuts (**on** blocks), and the **section** blocks: pieces
//! of the records lying on it. A trim residue lies on its owner's face
//! only (`IntersectionRecord::residue_owner_a`).

use std::collections::{BTreeMap, BTreeSet};

use popkern_core::topology::{BRep, FaceHandle};
use slotmap::Key;

use super::candidates::{CandidateVertex, Incidence};
use super::pave_blocks::{BlockOwner, PaveBlock};
use super::resolve::Resolution;
use crate::intersect::IntersectionRecord;

/// What the resolution put on one operand face.
#[derive(Debug, Clone, Default, PartialEq)]
pub(crate) struct FaceInfo {
    /// Vertices on the face's boundary, by group.
    pub on_vertices: BTreeSet<usize>,
    /// Vertices inside the face, by group.
    pub in_vertices: BTreeSet<usize>,
    /// Blocks of the face's own edges, by index into the blocks.
    pub on_blocks: Vec<usize>,
    /// Blocks of the records lying on the face, by index.
    pub section_blocks: Vec<usize>,
}

/// The edge and vertex handle bits on `face`'s boundary.
fn boundary(brep: &BRep, face: FaceHandle) -> (BTreeSet<u64>, BTreeSet<u64>) {
    let (mut edges, mut vertices) = (BTreeSet::new(), BTreeSet::new());
    let Some(f) = brep.faces.get(face) else {
        return (edges, vertices);
    };
    for lh in std::iter::once(f.outer_loop).chain(f.inner_loops.iter().copied()) {
        let Some(start) = brep.loops.get(lh).map(|l| l.half_edge) else {
            continue;
        };
        let mut he = start;
        for _ in 0..=brep.half_edges.len() {
            let Some(h) = brep.half_edges.get(he) else {
                break;
            };
            edges.insert(h.edge.data().as_ffi());
            vertices.insert(h.origin.data().as_ffi());
            he = h.next;
            if he == start {
                break;
            }
        }
    }
    (edges, vertices)
}

/// The faces a record lies on: (operand, face).
fn record_faces(r: &IntersectionRecord) -> Vec<(usize, FaceHandle)> {
    match r.residue_owner_a {
        Some(true) => vec![(0, r.face_a)],
        Some(false) => vec![(1, r.face_b)],
        None => vec![(0, r.face_a), (1, r.face_b)],
    }
}

/// FaceInfo for every operand face a record or a candidate lies on,
/// keyed by (operand, face).
pub(crate) fn face_infos(
    cands: &[CandidateVertex],
    res: &Resolution,
    operands: [&BRep; 2],
    records: &[IntersectionRecord],
    blocks: &[PaveBlock],
) -> BTreeMap<(usize, u64), FaceInfo> {
    let mut faces: BTreeMap<(usize, u64), FaceHandle> = BTreeMap::new();
    for r in records {
        for (o, f) in record_faces(r) {
            faces.insert((o, f.data().as_ffi()), f);
        }
    }
    for c in cands {
        for inc in &c.incidence {
            if let Incidence::Face { operand, face, .. } = *inc {
                faces.insert((operand, face.data().as_ffi()), face);
            }
        }
    }
    let mut out = BTreeMap::new();
    for (&(o, bits), &face) in &faces {
        let (edges, vertices) = boundary(operands[o], face);
        let mut info = FaceInfo::default();
        let mut inside = BTreeSet::new();
        for (i, c) in cands.iter().enumerate() {
            let group = res.group_of[i];
            for inc in &c.incidence {
                match *inc {
                    Incidence::Vertex { operand, vertex }
                        if operand == o && vertices.contains(&vertex.data().as_ffi()) =>
                    {
                        info.on_vertices.insert(group);
                    }
                    Incidence::Edge { operand, edge, .. }
                        if operand == o && edges.contains(&edge.data().as_ffi()) =>
                    {
                        info.on_vertices.insert(group);
                    }
                    Incidence::Face {
                        operand, face: f, ..
                    } if operand == o && f.data().as_ffi() == bits => {
                        inside.insert(group);
                    }
                    _ => {}
                }
            }
        }
        info.in_vertices = inside.difference(&info.on_vertices).copied().collect();
        for (k, b) in blocks.iter().enumerate() {
            match b.owner {
                BlockOwner::Edge(eo, e) if eo == o && edges.contains(&e) => info.on_blocks.push(k),
                BlockOwner::Record(r)
                    if records.get(r).is_some_and(|r| {
                        record_faces(r)
                            .iter()
                            .any(|&(ro, rf)| ro == o && rf.data().as_ffi() == bits)
                    }) =>
                {
                    info.section_blocks.push(k)
                }
                _ => {}
            }
        }
        out.insert((o, bits), info);
    }
    out
}
