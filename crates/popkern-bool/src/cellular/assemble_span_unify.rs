//! Coincident-span registry for the cellular assembler: two graph
//! pairs that connect the same two points with the same curve
//! midpoint are ONE brep edge (an intersection record that runs along
//! an edge the operand already has, or a boundary piece and its
//! overlap record on a flush interface).
//!
//! Coincidence is judged at the boolean's tolerance, not by node
//! identity: at the 1e-3 tolerance a marcher's record along a torus's
//! equator quarter starts at its own node 1.7e-5 from the original's
//! vertex, and an exact-node key unified four of eight quarters —
//! the annulus then consumed halves the bottom faces had claimed,
//! every ring re-wired reversed, and the served body kept two torus
//! faces with curveless edges (PK-420). An exact key is still tried
//! first (the tuned fast path); a tolerance scan over the registered
//! spans decides the rest.
//!
//! The scan compares GEOMETRY, not parameterisations (PK-777): a rod
//! cap's rim piece is an angle-parameterised arc while the flush box's
//! overlap record on the same arc is a rational quarter-circle span,
//! and their parameter midpoints coincide only on spans symmetric about
//! the quarter grid. With the rod's seam inside the overlap the two
//! midpoints sat 0.14 apart on one arc, the pairs became two edges each
//! used once, and the union declined non-manifold at every seam angle
//! off the axis line. Two spans with matching ends are one edge when
//! each one's midpoint lies on the other's piece; two distinct curves
//! through one node pair (a lens) fail that in both directions. The
//! on-piece band is the LINEAR tolerance, not the 10x endpoint gate: the
//! seam-straddling spans agree to 4e-15, while a cap plane 1e-6 past a
//! ring's extreme leaves a lens 1e-6 wide whose two curves sat inside
//! the 10x gate, were unified, and the union served 3.7 (1.3e-5) off
//! its volume (measured, tangent_boss_union at 8c35140c, against the
//! fuse volumes OCC computes 1e-3 and 2e-3 away; OCC's own fuse is
//! wrong that close). Since PK-780's closed-record projection that
//! lens no longer reaches the registry, so no fixture separates the
//! two bands today; the linear band stays as the coincidence rule.

use std::collections::HashMap;

use popkern_core::geometry::Point3;

use crate::merge::graph::IntersectionGraph;
use crate::merge::graph_duplicates::closest_on_piece;

/// Registered spans: `(start, end, midpoint, edge key)`.
pub(super) struct SpanRegistry {
    exact: HashMap<(usize, usize, i64, i64, i64), usize>,
    spans: Vec<(Point3, Point3, Point3, usize)>,
    /// Span indices by the grid cell (side 4·`gate`) of their first
    /// end: a span matching a query has that end within `gate` of one
    /// of the query's ends, so it sits in a cell that end's `gate` box
    /// overlaps. The scan over every span was quadratic in edge count
    /// (PK-802).
    cells: HashMap<(i64, i64, i64), Vec<usize>>,
    gate: f64,
    coincide: f64,
}

impl SpanRegistry {
    /// `gate` is the distance within which two spans' ends and
    /// midpoints are the same point; `coincide` (the linear tolerance)
    /// is the band within which one span's midpoint lies on the other's
    /// curve piece.
    pub(super) fn new(gate: f64, coincide: f64) -> Self {
        Self {
            exact: HashMap::new(),
            spans: Vec::new(),
            cells: HashMap::new(),
            gate,
            coincide,
        }
    }

    fn exact_key(nodes: (usize, usize), mid: Point3) -> (usize, usize, i64, i64, i64) {
        let (lo, hi) = (nodes.0.min(nodes.1), nodes.0.max(nodes.1));
        let q = |v: f64| (v / 1.0e-6).round() as i64;
        (lo, hi, q(mid.x), q(mid.y), q(mid.z))
    }

    /// The key of a registered span coincident with graph pair `key`'s
    /// span `(a, b, mid)`, in either direction.
    pub(super) fn prior_for(
        &self,
        graph: &IntersectionGraph,
        key: usize,
        (nodes, a, b, mid): ((usize, usize), Point3, Point3, Point3),
    ) -> Option<usize> {
        if let Some(&prior) = self.exact.get(&Self::exact_key(nodes, mid)) {
            return Some(prior);
        }
        let (g, on) = (self.gate, self.coincide);
        let mut candidates: Vec<usize> = Vec::new();
        for end in [a, b] {
            // The cells `end ± gate` overlaps (1-2 per axis).
            let lo = self.cell(&Point3::new(end.x - g, end.y - g, end.z - g));
            let hi = self.cell(&Point3::new(end.x + g, end.y + g, end.z + g));
            for x in lo.0..=hi.0 {
                for y in lo.1..=hi.1 {
                    for z in lo.2..=hi.2 {
                        if let Some(ids) = self.cells.get(&(x, y, z)) {
                            candidates.extend_from_slice(ids);
                        }
                    }
                }
            }
        }
        // Registration order, as the full scan visited them.
        candidates.sort_unstable();
        candidates.dedup();
        candidates
            .into_iter()
            .map(|i| &self.spans[i])
            .find(|(pa, pb, pm, prior)| {
                let same_ends = (pa.distance_to(&a) <= g && pb.distance_to(&b) <= g)
                    || (pa.distance_to(&b) <= g && pb.distance_to(&a) <= g);
                same_ends
                    && (pm.distance_to(&mid) <= g
                        || (piece_distance(graph, *prior, &mid) <= on
                            && piece_distance(graph, key, pm) <= on))
            })
            .map(|(_, _, _, prior)| *prior)
    }

    /// Register a span under `key`.
    pub(super) fn register(
        &mut self,
        nodes: (usize, usize),
        a: Point3,
        b: Point3,
        mid: Point3,
        key: usize,
    ) {
        self.exact.insert(Self::exact_key(nodes, mid), key);
        let cell = self.cell(&a);
        self.cells.entry(cell).or_default().push(self.spans.len());
        self.spans.push((a, b, mid, key));
    }

    /// The grid cell (side 4·`gate`) holding `p` (a non-positive gate
    /// collapses every point into one cell — the full scan).
    fn cell(&self, p: &Point3) -> (i64, i64, i64) {
        if self.gate.is_nan() || self.gate <= 0.0 {
            return (0, 0, 0);
        }
        let side = 4.0 * self.gate;
        let q = |v: f64| (v / side).floor() as i64;
        (q(p.x), q(p.y), q(p.z))
    }
}

/// The distance from `p` to graph pair `key`'s curve PIECE (its own
/// parameter window); a curveless pair is its straight chord.
fn piece_distance(graph: &IntersectionGraph, key: usize, p: &Point3) -> f64 {
    let he = &graph.half_edges[key];
    match he
        .curve
        .as_ref()
        .or_else(|| graph.half_edges[he.twin].curve.as_ref())
    {
        Some(curve) => closest_on_piece(curve, p).1,
        None => {
            let a = graph.nodes[he.origin].position;
            let b = graph.nodes[graph.half_edges[he.twin].origin].position;
            let ab = b - a;
            let len2 = ab.dot(&ab);
            let s = if len2 > 0.0 {
                ((*p - a).dot(&ab) / len2).clamp(0.0, 1.0)
            } else {
                0.0
            };
            (a + ab * s).distance_to(p)
        }
    }
}
