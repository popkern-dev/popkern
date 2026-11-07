//! Fin sense and the wedges about an edge (PK-1103).
//!
//! A fin is *forward* when it runs its edge the way the edge's first
//! fin does (`edge.half_edge`, which runs the stored curve forward).
//! The ring goes counter-clockwise about that direction seen from the
//! edge's end, so the space about an edge is cut by its faces into
//! wedges, one between each fin and its ring successor; a fin alone on
//! its edge (a free edge) bounds one wedge reaching round from its
//! face's front to its back.
//!
//! The side of a face facing a wedge follows from the fin's sense: a
//! face's effective normal points from its back to its front, its loop
//! runs counter-clockwise about that normal, so the wedge after a
//! forward fin lies on its face's front and the wedge before it on its
//! back (a reversed fin swaps them). Both sides of every wedge must be
//! in one shell (`validate_wedges_all`): that is what makes the face
//! sides of a shell a closed surface about its region, and it holds by
//! topology alone for a ring of two, whose two wedges are the same
//! whichever way the ring runs.

use super::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, ShellHandle};

/// One side of a face.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum FaceSide {
    /// The side the effective normal points away from.
    Back,
    /// The side the effective normal points into.
    Front,
}

impl FaceSide {
    /// The other side.
    pub fn opposite(self) -> Self {
        match self {
            Self::Back => Self::Front,
            Self::Front => Self::Back,
        }
    }
}

/// The space between two radially adjacent face uses of an edge,
/// bounded by `from`'s face on side `from_side` and `to`'s face on
/// side `to_side`. For a free edge `from == to` and the sides differ.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Wedge {
    /// The fin the wedge starts after (in ring order).
    pub from: HalfEdgeHandle,
    /// The side of `from`'s face facing the wedge.
    pub from_side: FaceSide,
    /// The fin the wedge ends at (`from`'s ring successor).
    pub to: HalfEdgeHandle,
    /// The side of `to`'s face facing the wedge.
    pub to_side: FaceSide,
}

impl BRep {
    /// Whether `fin` runs its edge the way the edge's first fin does.
    /// On an open edge that is reading its origin; on a closed edge
    /// (both ends one vertex) the first fin is forward and its mate in
    /// a ring of two is not, and any other fin of a longer ring is
    /// undecidable by topology (`None`), as is an unknown fin or one
    /// whose origin is neither end of its edge.
    pub fn fin_is_forward(&self, fin: HalfEdgeHandle) -> Option<bool> {
        let data = self.half_edges.get(fin)?;
        let first = self.edges.get(data.edge)?.half_edge;
        if fin == first {
            return Some(true);
        }
        let start = self.half_edges.get(first)?.origin;
        let end = self.fin_end(first)?;
        if start != end {
            return match data.origin {
                o if o == start => Some(true),
                o if o == end => Some(false),
                _ => None,
            };
        }
        (self.mate(first) == Some(fin)).then_some(false)
    }

    /// The wedges about `edge` in ring order from its first fin, or
    /// `None` when the edge has no face use (a wire edge, an unknown
    /// edge) or a fin's sense is undecidable ([`BRep::fin_is_forward`]).
    pub fn edge_wedges(&self, edge: EdgeHandle) -> Option<Vec<Wedge>> {
        let ring = self.radial_fins(edge);
        if ring.is_empty() {
            return None;
        }
        let senses: Vec<bool> = ring
            .iter()
            .map(|&f| self.fin_is_forward(f))
            .collect::<Option<_>>()?;
        let after = |forward: bool| {
            if forward {
                FaceSide::Front
            } else {
                FaceSide::Back
            }
        };
        let n = ring.len();
        Some(
            (0..n)
                .map(|i| {
                    let j = (i + 1) % n;
                    Wedge {
                        from: ring[i],
                        from_side: after(senses[i]),
                        to: ring[j],
                        to_side: after(senses[j]).opposite(),
                    }
                })
                .collect(),
        )
    }

    /// The face of fin `fin` (through its loop).
    pub fn fin_face(&self, fin: HalfEdgeHandle) -> Option<FaceHandle> {
        let loop_ = self.half_edges.get(fin)?.loop_;
        self.loops.get(loop_).map(|l| l.face)
    }

    /// The shell holding `face`'s `side`.
    pub fn face_side_shell(&self, face: FaceHandle, side: FaceSide) -> Option<ShellHandle> {
        self.faces.get(face).map(|f| match side {
            FaceSide::Back => f.back,
            FaceSide::Front => f.front,
        })
    }

    /// The shells holding a wedge's two bounding sides.
    pub fn wedge_shells(&self, wedge: &Wedge) -> Option<(ShellHandle, ShellHandle)> {
        let from = self.face_side_shell(self.fin_face(wedge.from)?, wedge.from_side)?;
        let to = self.face_side_shell(self.fin_face(wedge.to)?, wedge.to_side)?;
        Some((from, to))
    }
}
