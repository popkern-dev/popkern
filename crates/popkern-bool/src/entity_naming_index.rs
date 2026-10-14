//! Point index for the boolean's name matching (PK-803).
//!
//! Every output edge and vertex looked for its nearest input within a
//! `band` by scanning every input — quadratic in body size, ~4 ms of a
//! 64-slot disc's 64th cut. Inputs here are bucketed in a `band`-sized
//! grid: a match within `band` sits in one of the query's 27 cells.
//! Candidates come back in input order, so the first-minimum choice of
//! the full scan is kept.

use std::collections::HashMap;

use popkern_core::geometry::Point3;

/// Input points bucketed by `band`-sized cells.
pub(crate) struct NearIndex {
    cells: HashMap<(i64, i64, i64), Vec<usize>>,
    band: f64,
    count: usize,
}

impl NearIndex {
    /// Index `points` (their positions in the slice are the indices
    /// returned by [`NearIndex::near`]).
    pub(crate) fn new(points: impl Iterator<Item = Point3>, band: f64) -> Self {
        let mut index = Self {
            cells: HashMap::new(),
            band,
            count: 0,
        };
        for (i, p) in points.enumerate() {
            let cell = index.cell(&p);
            index.cells.entry(cell).or_default().push(i);
            index.count = i + 1;
        }
        index
    }

    fn cell(&self, p: &Point3) -> (i64, i64, i64) {
        if self.band.is_nan() || self.band <= 0.0 {
            return (0, 0, 0);
        }
        let q = |v: f64| (v / self.band).floor() as i64;
        (q(p.x), q(p.y), q(p.z))
    }

    /// Indices of every input that can lie within `band` of `p`, in
    /// input order (a superset; callers keep their own distance test).
    pub(crate) fn near(&self, p: &Point3) -> Vec<usize> {
        if !(p.x.is_finite() && p.y.is_finite() && p.z.is_finite()) {
            return (0..self.count).collect();
        }
        let (cx, cy, cz) = self.cell(p);
        let mut out = Vec::new();
        for dx in -1..=1_i64 {
            for dy in -1..=1_i64 {
                for dz in -1..=1_i64 {
                    let key = (
                        cx.saturating_add(dx),
                        cy.saturating_add(dy),
                        cz.saturating_add(dz),
                    );
                    if let Some(ids) = self.cells.get(&key) {
                        out.extend_from_slice(ids);
                    }
                }
            }
        }
        out.sort_unstable();
        out
    }
}
