//! Cell-structure checks of a general body (PK-1103, §5 checks 2, 4-6
//! and 8): every fin of an edge is on its ring, both sides of every
//! wedge lie in one shell, each shell is connected and its face sides
//! a closed surface, and the regions are well formed.
//!
//! The wedge check reads ring order and fin sense only (see
//! `topology::brep_wedges`): for a ring of two it is exact, for a
//! longer ring it trusts the stored order, whose agreement with the
//! geometry about the edge is checked apart. An edge with a fin of
//! undecidable sense (a closed edge of more than two fins) is not
//! judged, and neither are the shells reaching it.

use std::collections::{BTreeSet, HashMap};

use slotmap::Key;

use crate::error::{RegionDefect, ValidationError};
use crate::topology::{
    BRep, EdgeHandle, FaceHandle, FaceSide, HalfEdgeHandle, RegionKind, ShellHandle, VertexHandle,
};

/// A small union-find over dense indices.
struct Joins(Vec<usize>);

impl Joins {
    fn new(n: usize) -> Self {
        Self((0..n).collect())
    }

    fn root(&mut self, mut i: usize) -> usize {
        while self.0[i] != i {
            self.0[i] = self.0[self.0[i]];
            i = self.0[i];
        }
        i
    }

    fn join(&mut self, a: usize, b: usize) {
        let (a, b) = (self.root(a), self.root(b));
        self.0[a] = b;
    }
}

/// A side of a face, indexed densely: face index * 2 + side.
fn side_index(faces: &HashMap<FaceHandle, usize>, face: FaceHandle, side: FaceSide) -> usize {
    faces[&face] * 2 + usize::from(side == FaceSide::Front)
}

impl BRep {
    /// Check 2 (membership): every fin naming an edge is on the ring
    /// walked from the edge's first fin, and that first fin names the
    /// edge.
    pub fn validate_ring_membership_all(&self) -> Vec<ValidationError> {
        let mut errors = Vec::new();
        let mut on_ring: HashMap<EdgeHandle, BTreeSet<HalfEdgeHandle>> = HashMap::new();
        for (edge, data) in &self.edges {
            let first = data.half_edge;
            if !self.half_edges.get(first).is_some_and(|f| f.edge == edge) {
                errors.push(ValidationError::RingMembership {
                    edge,
                    half_edge: first,
                });
                continue;
            }
            let mut ring = BTreeSet::from([first]);
            let mut cursor = first;
            for _ in 0..self.half_edges.len() {
                match self.radial_next(cursor) {
                    Some(next) if ring.insert(next) => cursor = next,
                    _ => break,
                }
            }
            on_ring.insert(edge, ring);
        }
        for (fin, data) in &self.half_edges {
            if on_ring.get(&data.edge).is_some_and(|r| !r.contains(&fin)) {
                errors.push(ValidationError::RingMembership {
                    edge: data.edge,
                    half_edge: fin,
                });
            }
        }
        errors
    }

    /// Check 4: both sides bounding every wedge about every edge lie
    /// in one shell.
    pub fn validate_wedges_all(&self) -> Vec<ValidationError> {
        let mut errors = Vec::new();
        for edge in self.edges.keys() {
            for wedge in self.edge_wedges(edge).unwrap_or_default() {
                let Some((from_shell, to_shell)) = self.wedge_shells(&wedge) else {
                    continue; // reference checks own a missing loop or face
                };
                if from_shell != to_shell {
                    errors.push(ValidationError::WedgeShellMismatch {
                        edge,
                        from: wedge.from,
                        to: wedge.to,
                        from_shell,
                        to_shell,
                    });
                }
            }
        }
        errors
    }

    /// Checks 5 and 8: each shell's face sides, wire edges and acorn
    /// are one connected piece (sides through the wedges between them,
    /// wire edges through shared vertices and the vertices of the
    /// shell's faces), and its face sides, counted with an edge once
    /// per wedge of the shell and a vertex once per fan of the shell,
    /// satisfy V - E + 2F - L = 2 - 2g with g >= 0. Shells reaching a
    /// wedge that is undecidable or spans two shells are not judged
    /// (the wedge check reports the latter).
    pub fn validate_shell_cells_all(&self) -> Vec<ValidationError> {
        let face_index: HashMap<FaceHandle, usize> =
            self.faces.keys().enumerate().map(|(i, f)| (f, i)).collect();
        let fin_index: HashMap<HalfEdgeHandle, usize> = self
            .half_edges
            .keys()
            .enumerate()
            .map(|(i, h)| (h, i))
            .collect();
        let corner = |fin: HalfEdgeHandle, side: FaceSide| {
            fin_index[&fin] * 2 + usize::from(side == FaceSide::Front)
        };
        let mut sides = Joins::new(face_index.len() * 2);
        let mut corners = Joins::new(fin_index.len() * 2);
        let mut wedges_in: HashMap<ShellHandle, i64> = HashMap::new();
        let mut unjudged: BTreeSet<ShellHandle> = BTreeSet::new();
        for edge in self.edges.keys() {
            let ring = self.radial_fins(edge);
            if ring.is_empty() {
                continue;
            }
            let Some(wedges) = self.edge_wedges(edge) else {
                for f in ring.iter().filter_map(|&f| self.fin_face(f)) {
                    unjudged.extend([self.faces[f].back, self.faces[f].front]);
                }
                continue;
            };
            for w in wedges {
                let (Some(fa), Some(fb), Some((sa, sb))) = (
                    self.fin_face(w.from),
                    self.fin_face(w.to),
                    self.wedge_shells(&w),
                ) else {
                    continue;
                };
                if sa != sb {
                    unjudged.extend([sa, sb]);
                    continue;
                }
                *wedges_in.entry(sa).or_default() += 1;
                sides.join(
                    side_index(&face_index, fa, w.from_side),
                    side_index(&face_index, fb, w.to_side),
                );
                // The wedge joins the corners of its two sides at each
                // end of the edge: at the start a forward fin's corner
                // is its own, a reversed fin's its loop successor's.
                let ends = |fin: HalfEdgeHandle| -> Option<[HalfEdgeHandle; 2]> {
                    let next = self.half_edges.get(fin)?.next;
                    fin_index.get(&next)?;
                    Some(if self.fin_is_forward(fin)? {
                        [fin, next]
                    } else {
                        [next, fin]
                    })
                };
                if let (Some(a), Some(b)) = (ends(w.from), ends(w.to)) {
                    for k in 0..2 {
                        corners.join(corner(a[k], w.from_side), corner(b[k], w.to_side));
                    }
                }
            }
        }
        self.shell_findings(
            &face_index,
            &fin_index,
            &mut sides,
            &mut corners,
            &wedges_in,
            &unjudged,
        )
    }

    fn shell_findings(
        &self,
        face_index: &HashMap<FaceHandle, usize>,
        fin_index: &HashMap<HalfEdgeHandle, usize>,
        sides: &mut Joins,
        corners: &mut Joins,
        wedges_in: &HashMap<ShellHandle, i64>,
        unjudged: &BTreeSet<ShellHandle>,
    ) -> Vec<ValidationError> {
        let mut errors = Vec::new();
        // Per shell: its sides' roots, its corners' roots, F, L and
        // the vertices of empty loops.
        #[derive(Default)]
        struct Tally {
            pieces: BTreeSet<usize>,
            fans: BTreeSet<usize>,
            faces: i64,
            loops: i64,
            lone_vertices: i64,
        }
        let mut tally: HashMap<ShellHandle, Tally> = HashMap::new();
        let mut face_vertices: HashMap<FaceHandle, BTreeSet<VertexHandle>> = HashMap::new();
        for (face, data) in &self.faces {
            for side in [FaceSide::Back, FaceSide::Front] {
                let shell = self.face_side_shell(face, side).unwrap_or_default();
                let t = tally.entry(shell).or_default();
                t.pieces
                    .insert(sides.root(side_index(face_index, face, side)));
                t.faces += 1;
                for l in std::iter::once(data.outer_loop).chain(data.inner_loops.iter().copied()) {
                    t.loops += 1;
                    let Some(lp) = self.loops.get(l) else {
                        continue;
                    };
                    let mut fin = lp.half_edge;
                    if !fin_index.contains_key(&fin) {
                        t.lone_vertices += 1;
                        continue;
                    }
                    for _ in 0..fin_index.len() {
                        let Some(f) = self.half_edges.get(fin) else {
                            break;
                        };
                        let i = fin_index[&fin] * 2 + usize::from(side == FaceSide::Front);
                        t.fans.insert(corners.root(i));
                        face_vertices.entry(face).or_default().insert(f.origin);
                        fin = f.next;
                        if fin == lp.half_edge {
                            break;
                        }
                    }
                }
            }
        }
        for (shell, data) in &self.shells {
            let t = tally.remove(&shell).unwrap_or_default();
            let mut pieces = t.pieces.len();
            pieces += self.wire_pieces(shell, &data.wire_edges, &face_vertices, pieces > 0);
            pieces += usize::from(data.acorn.is_some());
            if pieces > 1 {
                errors.push(ValidationError::ShellDisconnected { shell, pieces });
            }
            if t.faces == 0 || unjudged.contains(&shell) {
                continue;
            }
            let v = t.fans.len() as i64 + t.lone_vertices;
            let e = wedges_in.get(&shell).copied().unwrap_or(0);
            let euler = v - e + 2 * t.faces - t.loops;
            if euler > 2 || euler % 2 != 0 {
                errors.push(ValidationError::ShellEuler { shell, euler });
            }
        }
        errors
    }

    /// How many pieces `wires` add to a shell: wire edges join through
    /// shared vertices, and a set of them touching a vertex of one of
    /// the shell's faces joins the face sides' piece (`has_sides`).
    fn wire_pieces(
        &self,
        shell: ShellHandle,
        wires: &[EdgeHandle],
        face_vertices: &HashMap<FaceHandle, BTreeSet<VertexHandle>>,
        has_sides: bool,
    ) -> usize {
        if wires.is_empty() {
            return 0;
        }
        let on_faces: BTreeSet<VertexHandle> = self
            .faces
            .iter()
            .filter(|(_, f)| f.back == shell || f.front == shell)
            .filter_map(|(h, _)| face_vertices.get(&h))
            .flatten()
            .copied()
            .collect();
        let mut joins = Joins::new(wires.len());
        let mut at: HashMap<VertexHandle, usize> = HashMap::new();
        let mut touches = vec![false; wires.len()];
        for (i, &e) in wires.iter().enumerate() {
            let Some((a, b)) = self.edge_vertices(e) else {
                continue;
            };
            for v in [a, b] {
                touches[i] |= on_faces.contains(&v);
                match at.get(&v) {
                    Some(&j) => joins.join(i, j),
                    None => {
                        at.insert(v, i);
                    }
                }
            }
        }
        let mut free = BTreeSet::new();
        let mut joined = BTreeSet::new();
        for (i, &touch) in touches.iter().enumerate() {
            let r = joins.root(i);
            if has_sides && touch {
                joined.insert(r);
            } else {
                free.insert(r);
            }
        }
        free.difference(&joined).count()
    }

    /// Check 6: the outside is a live void region with no outer shell;
    /// every other region has a live outer shell of its own; every
    /// shell's region is live.
    pub fn validate_regions_all(&self) -> Vec<ValidationError> {
        let mut errors = Vec::new();
        let defect = |region, defect| ValidationError::RegionStructure { region, defect };
        match self.regions.get(self.outside) {
            Some(r) if r.kind == RegionKind::Void => {
                if !r.shell.is_null() {
                    errors.push(defect(self.outside, RegionDefect::OutsideBounded));
                }
            }
            _ => errors.push(defect(self.outside, RegionDefect::OutsideNotVoid)),
        }
        for (region, data) in &self.regions {
            if region == self.outside {
                continue;
            }
            if !self
                .shells
                .get(data.shell)
                .is_some_and(|s| s.region == region)
            {
                errors.push(defect(region, RegionDefect::NoOuterShell));
            }
        }
        let missing: BTreeSet<_> = self
            .shells
            .values()
            .map(|s| s.region)
            .filter(|r| !self.regions.contains_key(*r))
            .collect();
        errors.extend(
            missing
                .into_iter()
                .map(|r| defect(r, RegionDefect::MissingRegion)),
        );
        errors
    }
}

#[cfg(test)]
#[path = "cells_tests.rs"]
mod tests;
