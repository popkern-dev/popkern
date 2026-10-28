//! Entities of the merged cell complex and its batch builder.
//!
//! Orientation conventions (fixed for every stage):
//!
//! - A face's loops wind counter-clockwise about the face's FOLDED
//!   normal (`SurfaceOnFace::folded_normal` of its primary owner),
//!   which is how the chart walk delivers fragment loops. The
//!   interior of the face lies to the LEFT of a fin seen from that
//!   normal: `folded_normal × fin_tangent` points into the face.
//! - `folded_is_outward` says whether that folded normal points out of
//!   the primary owner's material (`MaterialSense`); the face's
//!   [`SideIndex::Ahead`] side is the one its owner's OUTWARD normal
//!   points into.
//! - A fin is `forward` when it traverses its edge from `start` to
//!   `end`.

use popkern_core::geometry::{CurveOnEdge, Point3};
use popkern_core::topology::FaceHandle;

use super::ids::{EdgeId, FaceId, FinId, LoopId, ShellId, VertexId};
use super::operand_set::OperandSet;

/// A 0-cell.
#[derive(Debug, Clone)]
pub(crate) struct ComplexVertex {
    /// Position.
    pub position: Point3,
    /// Tolerance radius the vertex was committed at.
    pub tolerance: f64,
    /// The intersection-graph node this vertex lifts, if any.
    pub source_node: Option<usize>,
}

/// A 1-cell: any number of fins (2 for a manifold edge, 2m for a
/// non-manifold one, 0 for a wire edge).
#[derive(Debug, Clone)]
pub(crate) struct ComplexEdge {
    /// Start vertex.
    pub start: VertexId,
    /// End vertex (equal to `start` for a closed edge).
    pub end: VertexId,
    /// Curve, when the edge has one (degenerate pole edges may not).
    pub curve: Option<CurveOnEdge>,
    /// Tolerance of the edge's curve.
    pub tolerance: f64,
    /// Fins using this edge; angularly ordered once the radial stage
    /// has run (`radial_ordered`), arbitrary before.
    pub radial: Vec<FinId>,
    /// Whether `radial` is in angular order about the edge tangent.
    pub radial_ordered: bool,
    /// The intersection-graph edge key (`min(half-edge, twin)`) this
    /// edge lifts, if any.
    pub source_key: Option<usize>,
}

/// One face's use of one edge, inside one loop of that face (Lee &
/// Lee's partial edge). A seam edge carries two fins of the same face.
#[derive(Debug, Clone)]
pub(crate) struct Fin {
    /// The edge used.
    pub edge: EdgeId,
    /// The face using it.
    pub face: FaceId,
    /// The loop of `face` this fin lies in.
    pub loop_: LoopId,
    /// Traverses `edge` from `start` to `end`.
    pub forward: bool,
    /// Next fin in the loop.
    pub next: FinId,
    /// Previous fin in the loop.
    pub prev: FinId,
    /// The intersection-graph half-edge this fin lifts, if any.
    pub source_half_edge: Option<usize>,
}

/// A boundary loop of a face.
#[derive(Debug, Clone)]
pub(crate) struct ComplexLoop {
    /// Owning face.
    pub face: FaceId,
    /// A fin of the loop.
    pub first: FinId,
    /// Outer boundary (true) or hole.
    pub outer: bool,
}

/// One input face whose material the complex face carries.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct FaceOwner {
    /// Operand index.
    pub operand: usize,
    /// The input face.
    pub face: FaceHandle,
    /// Whether this owner's outward normal agrees with the complex
    /// face's (its primary owner's) outward normal. A same-domain face
    /// shared by a touching pair of opposite contact faces carries one
    /// owner of each sense.
    pub same_sense: bool,
}

/// Which side of a face, relative to its primary owner's outward
/// normal.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) enum SideIndex {
    /// The side the outward normal points into.
    Ahead = 0,
    /// The side behind the face (inside the primary owner).
    Behind = 1,
}

/// One side of a face (Lee & Lee's partial face).
#[derive(Debug, Clone, Default)]
pub(crate) struct FaceSide {
    /// The shell this side belongs to, once shells are assembled.
    pub shell: Option<ShellId>,
}

/// A 2-cell.
#[derive(Debug, Clone)]
pub(crate) struct ComplexFace {
    /// Input faces carried; the first is the primary owner, whose
    /// surface the face lies on.
    pub owners: Vec<FaceOwner>,
    /// Loops, outer first.
    pub loops: Vec<LoopId>,
    /// `[Ahead, Behind]`.
    pub sides: [FaceSide; 2],
    /// Whether the primary owner's folded normal points out of its
    /// material.
    pub folded_is_outward: bool,
    /// The fragments (walk outputs) this face lifts.
    pub source_fragments: Vec<usize>,
}

/// A connected set of face sides: one boundary component of a region.
#[derive(Debug, Clone, Default)]
pub(crate) struct ComplexShell {
    /// Its face sides, in discovery order.
    pub sides: Vec<(FaceId, SideIndex)>,
    /// Operands containing the region this shell bounds, once
    /// membership has been assigned.
    pub membership: Option<OperandSet>,
}

/// The merged cell complex.
#[derive(Debug, Clone, Default)]
pub(crate) struct CellComplex {
    /// 0-cells.
    pub vertices: Vec<ComplexVertex>,
    /// 1-cells.
    pub edges: Vec<ComplexEdge>,
    /// Face uses of edges.
    pub fins: Vec<Fin>,
    /// Face loops.
    pub loops: Vec<ComplexLoop>,
    /// 2-cells.
    pub faces: Vec<ComplexFace>,
    /// Region boundary components (empty until shells are assembled).
    pub shells: Vec<ComplexShell>,
    /// Number of operands imprinted.
    pub operand_count: usize,
}

impl CellComplex {
    /// An empty complex over `operand_count` operands.
    pub(crate) fn new(operand_count: usize) -> Self {
        Self {
            operand_count,
            ..Self::default()
        }
    }

    /// Add a vertex.
    pub(crate) fn add_vertex(
        &mut self,
        position: Point3,
        tolerance: f64,
        source_node: Option<usize>,
    ) -> VertexId {
        self.vertices.push(ComplexVertex {
            position,
            tolerance,
            source_node,
        });
        VertexId::from_index(self.vertices.len() - 1)
    }

    /// Add an edge with no fins yet.
    pub(crate) fn add_edge(
        &mut self,
        start: VertexId,
        end: VertexId,
        curve: Option<CurveOnEdge>,
        tolerance: f64,
        source_key: Option<usize>,
    ) -> EdgeId {
        self.edges.push(ComplexEdge {
            start,
            end,
            curve,
            tolerance,
            radial: Vec::new(),
            radial_ordered: false,
            source_key,
        });
        EdgeId::from_index(self.edges.len() - 1)
    }

    /// Add a face with no loops yet.
    pub(crate) fn add_face(
        &mut self,
        owners: Vec<FaceOwner>,
        folded_is_outward: bool,
        source_fragments: Vec<usize>,
    ) -> FaceId {
        self.faces.push(ComplexFace {
            owners,
            loops: Vec::new(),
            sides: [FaceSide::default(), FaceSide::default()],
            folded_is_outward,
            source_fragments,
        });
        FaceId::from_index(self.faces.len() - 1)
    }

    /// Add a loop to `face` from its fins in traversal order: each
    /// entry is `(edge, forward, source half-edge)`. Fins are linked
    /// cyclically and registered on their edges. An empty list adds
    /// nothing and returns `None`.
    pub(crate) fn add_loop(
        &mut self,
        face: FaceId,
        uses: &[(EdgeId, bool, Option<usize>)],
        outer: bool,
    ) -> Option<LoopId> {
        if uses.is_empty() {
            return None;
        }
        let loop_id = LoopId::from_index(self.loops.len());
        let base = self.fins.len();
        let n = uses.len();
        for (i, &(edge, forward, source_half_edge)) in uses.iter().enumerate() {
            let fin = FinId::from_index(base + i);
            self.fins.push(Fin {
                edge,
                face,
                loop_: loop_id,
                forward,
                next: FinId::from_index(base + (i + 1) % n),
                prev: FinId::from_index(base + (i + n - 1) % n),
                source_half_edge,
            });
            if let Some(e) = self.edges.get_mut(edge.index()) {
                e.radial.push(fin);
                e.radial_ordered = false;
            }
        }
        self.loops.push(ComplexLoop {
            face,
            first: FinId::from_index(base),
            outer,
        });
        if let Some(f) = self.faces.get_mut(face.index()) {
            if outer {
                f.loops.insert(0, loop_id);
            } else {
                f.loops.push(loop_id);
            }
        }
        Some(loop_id)
    }

    /// The vertex a fin leaves.
    pub(crate) fn fin_start(&self, fin: FinId) -> Option<VertexId> {
        let f = self.fins.get(fin.index())?;
        let e = self.edges.get(f.edge.index())?;
        Some(if f.forward { e.start } else { e.end })
    }

    /// The vertex a fin arrives at.
    pub(crate) fn fin_end(&self, fin: FinId) -> Option<VertexId> {
        let f = self.fins.get(fin.index())?;
        let e = self.edges.get(f.edge.index())?;
        Some(if f.forward { e.end } else { e.start })
    }

    /// The fins of a loop in traversal order (bounded by the fin count,
    /// so a corrupt cycle cannot loop forever).
    pub(crate) fn loop_fins(&self, loop_id: LoopId) -> Vec<FinId> {
        let Some(l) = self.loops.get(loop_id.index()) else {
            return Vec::new();
        };
        let mut out = vec![l.first];
        let mut cursor = l.first;
        for _ in 0..self.fins.len() {
            let Some(f) = self.fins.get(cursor.index()) else {
                break;
            };
            if f.next == l.first {
                break;
            }
            cursor = f.next;
            out.push(cursor);
        }
        out
    }
}
