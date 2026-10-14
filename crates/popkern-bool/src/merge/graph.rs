//! Intersection graph data structures for the graph-based merge algorithm.
//!
//! The intersection graph represents the combined topology of two BRep operands
//! plus their intersection curves, using indexed arrays of nodes and half-edges.

// These types are used only by sibling modules and tests for now;
// they will be consumed by later pipeline stages.
#![allow(dead_code)]

use std::collections::{BTreeMap, HashMap};

use popkern_core::{
    geometry::{CurveOnEdge, Point3},
    topology::{EdgeHandle, FaceHandle, VertexHandle},
};

/// Identifies which boolean operand a piece of topology comes from.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) enum Operand {
    /// The first (left-hand) operand.
    A,
    /// The second (right-hand) operand.
    B,
}

/// A face qualified by which operand it belongs to.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct FaceRef {
    /// Which operand this face belongs to.
    pub operand: Operand,
    /// Handle to the face in its source BRep.
    pub face: FaceHandle,
}

/// Tracks the origin of a graph node.
#[derive(Debug, Clone)]
pub(crate) enum NodeSource {
    /// Copied from an original vertex in operand A.
    OriginalA(VertexHandle),
    /// Copied from an original vertex in operand B.
    OriginalB(VertexHandle),
    /// Created at a curve-curve or edge-face intersection point.
    Intersection,
}

/// Tracks the origin of a graph edge.
#[derive(Debug, Clone)]
pub(crate) enum EdgeSource {
    /// Copied from an original edge in operand A.
    OriginalA(EdgeHandle),
    /// Copied from an original edge in operand B.
    OriginalB(EdgeHandle),
    /// A shared edge lying on the intersection of two faces.
    Shared {
        /// Face from operand A.
        face_a: FaceHandle,
        /// Face from operand B.
        face_b: FaceHandle,
    },
    /// Created by splitting an original edge from operand A.
    SplitA {
        /// The original edge that was split.
        original: EdgeHandle,
    },
    /// Created by splitting an original edge from operand B.
    SplitB {
        /// The original edge that was split.
        original: EdgeHandle,
    },
}

/// A node (vertex) in the intersection graph.
#[derive(Debug)]
pub(crate) struct GraphNode {
    /// 3D position.
    pub position: Point3,
    /// UV coordinates of this node on each adjacent face.
    pub uv_coords: HashMap<FaceRef, (f64, f64)>,
    /// Where this node came from.
    pub source: NodeSource,
    /// Indices of outgoing half-edges in the graph's `half_edges` array.
    pub outgoing: Vec<usize>,
}

/// A directed half-edge in the intersection graph.
#[derive(Debug)]
pub(crate) struct GraphHalfEdge {
    /// Index of the origin node.
    pub origin: usize,
    /// Index of the twin (opposite-direction) half-edge.
    pub twin: usize,
    /// The face this half-edge borders, if assigned.
    pub face: Option<FaceRef>,
    /// Geometric curve along this half-edge, if any.
    pub curve: Option<CurveOnEdge>,
    /// UV parameter-space curve (pcurve) for the face this half-edge borders.
    pub pcurve: Option<popkern_core::geometry::pcurve::Pcurve>,
    /// Where this edge came from.
    pub source: EdgeSource,
}

/// The intersection graph: a combined representation of two operands'
/// topology plus intersection curves, ready for loop walking and face
/// reconstruction.
#[derive(Debug)]
pub(crate) struct IntersectionGraph {
    /// All nodes (vertices) in the graph.
    pub nodes: Vec<GraphNode>,
    /// All directed half-edges in the graph.
    pub half_edges: Vec<GraphHalfEdge>,
    /// Maps each face to the half-edge indices that border it.
    /// `BTreeMap` (not `HashMap`): `walk_all_fragments` iterates this
    /// map, and fragment order determines result face-handle order —
    /// unordered iteration made boolean output face handles vary
    /// run-to-run (PK-099 determinism hunt).
    pub face_edges: BTreeMap<FaceRef, Vec<usize>>,
    /// Shared half-edges that took over an operand SEAM on a face: the
    /// duplicate-cover pass removed the seam's own half-edges in their
    /// favour, so the chart treats them as that face's seam copies
    /// (PK-750; see `graph_duplicates`).
    pub seam_substitutes: std::collections::BTreeSet<(FaceRef, usize)>,
    /// Node indices keyed by their x coordinate (order-preserving bit
    /// key), so `find_node` scans a slab instead of every node: graph
    /// build looked up every operand vertex against every node before
    /// it, quadratic in body size (PK-802). Positions are never moved
    /// after `add_node`.
    node_index: BTreeMap<u64, Vec<usize>>,
}

/// Order-preserving key of an `f64` (total order: −0 < +0, NaN last).
fn order_key(x: f64) -> u64 {
    let bits = x.to_bits();
    if bits >> 63 == 1 {
        !bits
    } else {
        bits | (1 << 63)
    }
}

impl IntersectionGraph {
    /// Create an empty intersection graph.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            half_edges: Vec::new(),
            face_edges: BTreeMap::new(),
            seam_substitutes: std::collections::BTreeSet::new(),
            node_index: BTreeMap::new(),
        }
    }

    /// Add a node at the given position with the given source.
    /// Returns the index of the new node.
    pub fn add_node(&mut self, position: Point3, source: NodeSource) -> usize {
        let idx = self.nodes.len();
        self.node_index
            .entry(order_key(position.x))
            .or_default()
            .push(idx);
        self.nodes.push(GraphNode {
            position,
            uv_coords: HashMap::new(),
            source,
            outgoing: Vec::new(),
        });
        idx
    }

    /// Add a half-edge pair (forward + twin) between two nodes.
    ///
    /// Both half-edges share the same `source` and `curve`. The twin's curve
    /// is set to `None` (the forward half-edge owns the geometry).
    /// Returns `(forward_idx, twin_idx)`.
    pub fn add_half_edge_pair(
        &mut self,
        start: usize,
        end: usize,
        source: EdgeSource,
        curve: Option<CurveOnEdge>,
    ) -> (usize, usize) {
        self.add_half_edge_pair_with_pcurves(start, end, source, curve, None, None)
    }

    /// Add a half-edge pair with explicit pcurves for each direction.
    ///
    /// Like [`add_half_edge_pair`](Self::add_half_edge_pair), but also
    /// attaches UV parameter-space curves to the forward and twin
    /// half-edges respectively.
    pub fn add_half_edge_pair_with_pcurves(
        &mut self,
        start: usize,
        end: usize,
        source: EdgeSource,
        curve: Option<CurveOnEdge>,
        pcurve_fwd: Option<popkern_core::geometry::pcurve::Pcurve>,
        pcurve_twin: Option<popkern_core::geometry::pcurve::Pcurve>,
    ) -> (usize, usize) {
        let fwd_idx = self.half_edges.len();
        let twin_idx = fwd_idx + 1;

        self.half_edges.push(GraphHalfEdge {
            origin: start,
            twin: twin_idx,
            face: None,
            curve,
            pcurve: pcurve_fwd,
            source: source.clone(),
        });
        self.half_edges.push(GraphHalfEdge {
            origin: end,
            twin: fwd_idx,
            face: None,
            curve: None,
            pcurve: pcurve_twin,
            source,
        });

        self.nodes[start].outgoing.push(fwd_idx);
        self.nodes[end].outgoing.push(twin_idx);

        (fwd_idx, twin_idx)
    }

    /// Find a node whose position is within `tolerance` of `position`.
    /// Returns `None` if no such node exists. O(N) scan.
    ///
    /// Uses [`popkern_core::adaptive_precision::adaptive_dist_sq_lt`]
    /// for the distance comparison — avoids the sqrt in
    /// [`Point3::distance_to`] and escalates to double-precision
    /// when the f64 answer falls inside the ambiguity margin, so
    /// snap decisions at the edge of the tolerance envelope are
    /// deterministic and order-independent.
    pub fn find_node(&self, position: &Point3, tolerance: f64) -> Option<usize> {
        let threshold_sq = tolerance * tolerance;
        let p = (position.x, position.y, position.z);
        // The slab is widened past the predicate's reach (relative
        // and absolute slack), so every node the exact distance test
        // accepts is in it; the lowest accepted index is the linear
        // scan's answer.
        let reach = tolerance.abs() * (1.0 + 1e-6) + position.x.abs() * 1e-12 + 1e-300;
        let (lo, hi) = (order_key(position.x - reach), order_key(position.x + reach));
        if !(lo <= hi) || !tolerance.is_finite() {
            return self.find_node_linear(p, threshold_sq);
        }
        self.node_index
            .range(lo..=hi)
            .flat_map(|(_, ids)| ids.iter().copied())
            .filter(|&i| {
                let n = &self.nodes[i];
                popkern_core::adaptive_precision::adaptive_dist_sq_lt(
                    (n.position.x, n.position.y, n.position.z),
                    p,
                    threshold_sq,
                )
            })
            .min()
    }

    /// `find_node` by scanning every node (non-finite queries).
    fn find_node_linear(&self, p: (f64, f64, f64), threshold_sq: f64) -> Option<usize> {
        self.nodes.iter().position(|n| {
            popkern_core::adaptive_precision::adaptive_dist_sq_lt(
                (n.position.x, n.position.y, n.position.z),
                p,
                threshold_sq,
            )
        })
    }

    /// Assign a face to a half-edge and update the `face_edges` index.
    pub fn set_half_edge_face(&mut self, he_idx: usize, face: FaceRef) {
        self.half_edges[he_idx].face = Some(face);
        self.face_edges.entry(face).or_default().push(he_idx);
    }
}

/// The `uv-scramble` trace topic: image-honesty sweep over the
/// graph's shared-record pcurves — every registered pcurve, evaluated
/// on its face's surface, must trace its own 3D curve (either
/// direction). A pcurve imaging off its curve at feature scale is the
/// mixed-branch / wrong-end class the mesher later tears on
/// (PK-485). Call between graph-mutation stages with a stage label.
pub(crate) fn trace_uv_scramble(
    stage: &str,
    graph: &IntersectionGraph,
    brep_a: &popkern_core::topology::BRep,
    brep_b: &popkern_core::topology::BRep,
    tolerance: &popkern_core::tolerance::Tolerance,
) {
    if !popkern_core::trace::trace_enabled("uv-scramble") {
        return;
    }
    use popkern_core::geometry::curve::CurveEval as _;
    use popkern_core::geometry::surface::SurfaceEval as _;
    for (i, he) in graph.half_edges.iter().enumerate() {
        let (Some(c), Some(pc), Some(fr)) = (he.curve.as_ref(), he.pcurve.as_ref(), &he.face)
        else {
            continue;
        };
        let brep = match fr.operand {
            Operand::A => brep_a,
            Operand::B => brep_b,
        };
        let Some(surf) = brep.faces.get(fr.face).and_then(|f| f.surface.as_ref()) else {
            continue;
        };
        let (p0, p1) = pc.parameter_range();
        let mut worst = 0.0_f64;
        let mut rev = 0.0_f64;
        for k in 0..=8 {
            let f = k as f64 / 8.0;
            let q = pc.point_at(p0 + (p1 - p0) * f);
            let ps = surf.geometry.point_at(q.x, q.y);
            let fwd3 = c.geometry.point_at(c.t_start + (c.t_end - c.t_start) * f);
            let rev3 = c.geometry.point_at(c.t_end + (c.t_start - c.t_end) * f);
            worst = worst.max(ps.distance_to(&fwd3));
            rev = rev.max(ps.distance_to(&rev3));
        }
        if worst.min(rev) > tolerance.linear() * 100.0 {
            let a = graph.nodes[he.origin].position;
            let b = graph.nodes[graph.half_edges[he.twin].origin].position;
            eprintln!(
                "  [uv-scramble] {stage}: graph he {i} {fr:?} src {:?} images {worst:.3e} fwd / {rev:.3e} rev ({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3})",
                he.source, a.x, a.y, a.z, b.x, b.y, b.z
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_graph() {
        let g = IntersectionGraph::new();
        assert!(g.nodes.is_empty());
        assert!(g.half_edges.is_empty());
        assert!(g.face_edges.is_empty());
    }

    #[test]
    fn add_node_returns_index() {
        let mut g = IntersectionGraph::new();
        let i0 = g.add_node(Point3::new(0.0, 0.0, 0.0), NodeSource::Intersection);
        let i1 = g.add_node(Point3::new(1.0, 0.0, 0.0), NodeSource::Intersection);
        assert_eq!(i0, 0);
        assert_eq!(i1, 1);
        assert_eq!(g.nodes.len(), 2);
    }

    #[test]
    fn add_half_edge_pair_wiring() {
        let mut g = IntersectionGraph::new();
        let n0 = g.add_node(Point3::new(0.0, 0.0, 0.0), NodeSource::Intersection);
        let n1 = g.add_node(Point3::new(1.0, 0.0, 0.0), NodeSource::Intersection);

        let edge_h = popkern_core::topology::EdgeHandle::default();
        let (fwd, twin) = g.add_half_edge_pair(n0, n1, EdgeSource::OriginalA(edge_h), None);

        // Twin pairing
        assert_eq!(g.half_edges[fwd].twin, twin);
        assert_eq!(g.half_edges[twin].twin, fwd);

        // Origins
        assert_eq!(g.half_edges[fwd].origin, n0);
        assert_eq!(g.half_edges[twin].origin, n1);

        // Outgoing lists
        assert!(g.nodes[n0].outgoing.contains(&fwd));
        assert!(g.nodes[n1].outgoing.contains(&twin));
    }

    /// The x-slab index answers exactly as the linear scan (lowest
    /// accepted index), including clustered, signed-zero and
    /// tolerance-edge queries.
    #[test]
    fn find_node_index_matches_linear_scan() {
        let mut g = IntersectionGraph::new();
        let mut seed = 0x2545_f491_4f6c_dd1du64;
        let mut next = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            (seed >> 11) as f64 / (1u64 << 53) as f64
        };
        for i in 0..400 {
            // Clusters within 1e-6 of each other and of ±0.
            let c = (i % 7) as f64 - 3.0;
            let p = Point3::new(c * 1e-7 * next(), c + next(), 2.0 * next());
            g.add_node(p, NodeSource::Intersection);
            if i % 5 == 0 {
                g.add_node(Point3::new(-0.0, p.y, p.z), NodeSource::Intersection);
            }
        }
        let positions: Vec<Point3> = g.nodes.iter().map(|n| n.position).collect();
        for (k, p) in positions.iter().enumerate() {
            for tol in [0.0, 1e-9, 1e-6, 1e-3, 0.3, 5.0] {
                let q = Point3::new(p.x + tol * 0.999 * (k % 3) as f64 - tol, p.y, p.z);
                let linear = g.find_node_linear((q.x, q.y, q.z), tol * tol);
                assert_eq!(g.find_node(&q, tol), linear, "query {q:?} tol {tol}");
            }
        }
    }

    #[test]
    fn find_node_within_tolerance() {
        let mut g = IntersectionGraph::new();
        g.add_node(Point3::new(1.0, 2.0, 3.0), NodeSource::Intersection);

        let found = g.find_node(&Point3::new(1.0, 2.0, 3.0 + 1e-8), 1e-6);
        assert_eq!(found, Some(0));

        let not_found = g.find_node(&Point3::new(10.0, 0.0, 0.0), 1e-6);
        assert!(not_found.is_none());
    }

    #[test]
    fn set_half_edge_face_updates_index() {
        let mut g = IntersectionGraph::new();
        let n0 = g.add_node(Point3::new(0.0, 0.0, 0.0), NodeSource::Intersection);
        let n1 = g.add_node(Point3::new(1.0, 0.0, 0.0), NodeSource::Intersection);
        let edge_h = popkern_core::topology::EdgeHandle::default();
        let (fwd, _twin) = g.add_half_edge_pair(n0, n1, EdgeSource::OriginalA(edge_h), None);

        let face_ref = FaceRef {
            operand: Operand::A,
            face: popkern_core::topology::FaceHandle::default(),
        };
        g.set_half_edge_face(fwd, face_ref);

        assert_eq!(g.half_edges[fwd].face, Some(face_ref));
        assert_eq!(g.face_edges[&face_ref], vec![fwd]);
    }
}
