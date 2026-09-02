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

use std::collections::HashMap;

use popkern_core::geometry::Point3;

/// Registered spans: `(start, end, midpoint, edge key)`.
pub(super) struct SpanRegistry {
    exact: HashMap<(usize, usize, i64, i64, i64), usize>,
    spans: Vec<(Point3, Point3, Point3, usize)>,
    gate: f64,
}

impl SpanRegistry {
    /// `gate` is the distance within which two spans' ends and
    /// midpoints are the same point.
    pub(super) fn new(gate: f64) -> Self {
        Self {
            exact: HashMap::new(),
            spans: Vec::new(),
            gate,
        }
    }

    fn exact_key(nodes: (usize, usize), mid: Point3) -> (usize, usize, i64, i64, i64) {
        let (lo, hi) = (nodes.0.min(nodes.1), nodes.0.max(nodes.1));
        let q = |v: f64| (v / 1.0e-6).round() as i64;
        (lo, hi, q(mid.x), q(mid.y), q(mid.z))
    }

    /// The key of a registered span coincident with `(a, b, mid)`, in
    /// either direction.
    pub(super) fn prior_for(
        &self,
        nodes: (usize, usize),
        a: Point3,
        b: Point3,
        mid: Point3,
    ) -> Option<usize> {
        if let Some(&key) = self.exact.get(&Self::exact_key(nodes, mid)) {
            return Some(key);
        }
        let g = self.gate;
        self.spans
            .iter()
            .find(|(pa, pb, pm, _)| {
                pm.distance_to(&mid) <= g
                    && ((pa.distance_to(&a) <= g && pb.distance_to(&b) <= g)
                        || (pa.distance_to(&b) <= g && pb.distance_to(&a) <= g))
            })
            .map(|(_, _, _, key)| *key)
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
        self.spans.push((a, b, mid, key));
    }
}
