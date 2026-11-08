//! Joining and separating the shells of a region as gluing connects
//! its boundary pieces or ungluing disconnects them.
//!
//! A shell is one connected piece of a region's boundary: face sides
//! joined through the wedges about their edges (both sides of a wedge
//! lie in one shell), with the wire edges reaching them. Gluing two
//! edges makes new wedges whose sides may lie in two shells of one
//! region, which then join; ungluing can leave a shell in pieces, each
//! then a shell of its own.

use std::collections::{BTreeMap, BTreeSet, HashMap};

use crate::topology::{BRep, FaceHandle, FaceSide, Shell, ShellHandle};

impl BRep {
    /// Move everything `other` holds into `target` (both of one
    /// region): face sides, wire edges, and the region's anchor; `other`
    /// goes.
    pub(super) fn merge_shell_into(&mut self, target: ShellHandle, other: ShellHandle) {
        if target == other || !self.shells.contains_key(other) {
            return;
        }
        for face in self.faces.values_mut() {
            if face.back == other {
                face.back = target;
            }
            if face.front == other {
                face.front = target;
            }
        }
        let wires = std::mem::take(&mut self.shells[other].wire_edges);
        self.shells[target].wire_edges.extend(wires);
        let region = self.shells[other].region;
        if self.regions.get(region).is_some_and(|r| r.shell == other) {
            self.regions[region].shell = target;
        }
        if !self.faces.contains_key(self.shells[target].face) {
            self.shells[target].face = self.shells[other].face;
        }
        self.shells.remove(other);
    }

    /// The pieces of `shell`'s face sides joined through wedges, each a
    /// set of (face, side), ordered by their least member. An edge whose
    /// wedges are undecidable joins every side of the shell on it.
    pub(super) fn shell_pieces(&self, shell: ShellHandle) -> Vec<BTreeSet<(FaceHandle, FaceSide)>> {
        let sides: Vec<(FaceHandle, FaceSide)> = self
            .faces
            .iter()
            .flat_map(|(h, f)| {
                let mut v = Vec::new();
                if f.back == shell {
                    v.push((h, FaceSide::Back));
                }
                if f.front == shell {
                    v.push((h, FaceSide::Front));
                }
                v
            })
            .collect();
        let index: HashMap<(FaceHandle, FaceSide), usize> =
            sides.iter().enumerate().map(|(i, &s)| (s, i)).collect();
        let mut parent: Vec<usize> = (0..sides.len()).collect();
        fn root(p: &mut [usize], mut i: usize) -> usize {
            while p[i] != i {
                p[i] = p[p[i]];
                i = p[i];
            }
            i
        }
        let join = |p: &mut Vec<usize>, a: usize, b: usize| {
            let (x, y) = (root(p, a), root(p, b));
            p[x] = y;
        };
        let faces: BTreeSet<FaceHandle> = sides.iter().map(|s| s.0).collect();
        let edges: BTreeSet<_> = self
            .half_edges
            .values()
            .filter(|h| {
                self.loops
                    .get(h.loop_)
                    .is_some_and(|l| faces.contains(&l.face))
            })
            .map(|h| h.edge)
            .collect();
        for edge in edges {
            match self.edge_wedges(edge) {
                Some(wedges) => {
                    for w in wedges {
                        let (Some(fa), Some(fb)) = (self.fin_face(w.from), self.fin_face(w.to))
                        else {
                            continue;
                        };
                        if let (Some(&a), Some(&b)) =
                            (index.get(&(fa, w.from_side)), index.get(&(fb, w.to_side)))
                        {
                            join(&mut parent, a, b);
                        }
                    }
                }
                None => {
                    let on: Vec<usize> = self
                        .radial_fins(edge)
                        .iter()
                        .filter_map(|&f| self.fin_face(f))
                        .flat_map(|f| [(f, FaceSide::Back), (f, FaceSide::Front)])
                        .filter_map(|s| index.get(&s).copied())
                        .collect();
                    for w in on.windows(2) {
                        join(&mut parent, w[0], w[1]);
                    }
                }
            }
        }
        let mut groups: BTreeMap<usize, BTreeSet<(FaceHandle, FaceSide)>> = BTreeMap::new();
        for (i, &s) in sides.iter().enumerate() {
            let r = root(&mut parent, i);
            groups.entry(r).or_default().insert(s);
        }
        let mut out: Vec<_> = groups.into_values().collect();
        out.sort_by(|a, b| a.first().cmp(&b.first()));
        out
    }

    /// Give every piece of `shell` but the one holding its anchor a
    /// shell of its own in the same region, wire edges going with the
    /// piece whose faces they reach (or staying). Returns the new
    /// shells.
    pub(super) fn separate_shell_pieces(&mut self, shell: ShellHandle) -> Vec<ShellHandle> {
        let pieces = self.shell_pieces(shell);
        if pieces.len() < 2 {
            return Vec::new();
        }
        let anchor = self.shells[shell].face;
        let region = self.shells[shell].region;
        let stay = pieces
            .iter()
            .position(|p| p.iter().any(|s| s.0 == anchor))
            .unwrap_or(0);
        let mut made = Vec::new();
        for (i, piece) in pieces.into_iter().enumerate() {
            if i == stay {
                continue;
            }
            let first = piece.first().map_or(FaceHandle::default(), |s| s.0);
            let new = self.shells.insert(Shell::new(first, region));
            for (face, side) in piece {
                let f = &mut self.faces[face];
                match side {
                    FaceSide::Back => f.back = new,
                    FaceSide::Front => f.front = new,
                }
            }
            made.push(new);
        }
        // Wire edges follow their pieces.
        for &new in &made {
            let vertices = self.shell_face_vertices(new);
            let wires = self.shells[shell].wire_edges.clone();
            let (moved, kept): (Vec<_>, Vec<_>) = wires.into_iter().partition(|&e| {
                self.edge_vertices(e)
                    .is_some_and(|(a, b)| vertices.contains(&a) || vertices.contains(&b))
            });
            self.shells[shell].wire_edges = kept;
            self.shells[new].wire_edges.extend(moved);
        }
        made
    }

    /// The vertices of the faces with a side in `shell`.
    fn shell_face_vertices(&self, shell: ShellHandle) -> BTreeSet<crate::topology::VertexHandle> {
        self.half_edges
            .values()
            .filter(|h| {
                self.loops
                    .get(h.loop_)
                    .and_then(|l| self.faces.get(l.face))
                    .is_some_and(|f| f.back == shell || f.front == shell)
            })
            .map(|h| h.origin)
            .collect()
    }
}
