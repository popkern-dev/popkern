//! Spatial index for the crack zipper's on-curve location (split from
//! `crack_zipper.rs`, 500-line cap).
//!
//! Location asks, for one mesh vertex, which two-sided edge polyline
//! passes within that edge's band and where along it. Scanning every
//! segment of every edge per vertex is quadratic in model size
//! (measured: turbine.stp, 131 cones, never finished meshing once its
//! geometry imported correctly and thousands of rim stations went to
//! the zipper). A segment can only accept `p` if `p` lies inside the
//! segment's bounding box grown by the edge's band, so hashing those
//! grown boxes onto a uniform grid returns exactly the candidates the
//! full scan could accept. Candidates are visited in (edge, segment)
//! order with a strict `<`, the full scan's tie-break.

use std::collections::HashMap;

use popkern_core::geometry::Point3;

use super::{EdgePolyline, project_point_segment};

/// Grid cells → the `(edge, segment)` pairs whose grown boxes touch them.
pub(super) struct SegmentIndex {
    cell: f64,
    cells: HashMap<(i64, i64, i64), Vec<(u32, u32)>>,
}

impl SegmentIndex {
    /// Index every segment of `edges`.
    pub(super) fn build(edges: &[EdgePolyline]) -> SegmentIndex {
        let mut lengths: Vec<f64> = edges
            .iter()
            .flat_map(|e| e.pts.windows(2).map(|w| w[0].distance_to(&w[1])))
            .collect();
        let max_band = edges.iter().map(|e| e.band).fold(0.0, f64::max);
        lengths.sort_by(f64::total_cmp);
        let median = lengths.get(lengths.len() / 2).copied().unwrap_or(0.0);
        let longest = lengths.last().copied().unwrap_or(0.0);
        // Cells no smaller than the band (a query touches few), the
        // typical segment, or a sixteenth of the longest one (a long
        // segment touches a bounded number of cells).
        let cell = (4.0 * max_band).max(median).max(longest / 16.0).max(1e-9);
        let mut cells: HashMap<(i64, i64, i64), Vec<(u32, u32)>> = HashMap::new();
        for (ei, edge) in edges.iter().enumerate() {
            for (si, w) in edge.pts.windows(2).enumerate() {
                let lo = key(
                    cell,
                    w[0].x.min(w[1].x) - edge.band,
                    w[0].y.min(w[1].y) - edge.band,
                    w[0].z.min(w[1].z) - edge.band,
                );
                let hi = key(
                    cell,
                    w[0].x.max(w[1].x) + edge.band,
                    w[0].y.max(w[1].y) + edge.band,
                    w[0].z.max(w[1].z) + edge.band,
                );
                for i in lo.0..=hi.0 {
                    for j in lo.1..=hi.1 {
                        for k in lo.2..=hi.2 {
                            cells
                                .entry((i, j, k))
                                .or_default()
                                .push((ei as u32, si as u32));
                        }
                    }
                }
            }
        }
        SegmentIndex { cell, cells }
    }

    /// Best on-curve location of `p` across all polylines within band:
    /// `(edge index, arc position)`.
    pub(super) fn locate(&self, edges: &[EdgePolyline], p: &Point3) -> Option<(usize, f64)> {
        let candidates = self.cells.get(&key(self.cell, p.x, p.y, p.z))?;
        let mut best: Option<(usize, f64, f64)> = None;
        for &(ei, si) in candidates {
            let (ei, si) = (ei as usize, si as usize);
            let edge = &edges[ei];
            let (f, foot) = project_point_segment(p, &edge.pts[si], &edge.pts[si + 1]);
            let d = foot.distance_to(p);
            if d <= edge.band && best.is_none_or(|(_, _, bd)| d < bd) {
                best = Some((ei, edge.cum[si] + (edge.cum[si + 1] - edge.cum[si]) * f, d));
            }
        }
        best.map(|(ei, s, _)| (ei, s))
    }
}

fn key(cell: f64, x: f64, y: f64, z: f64) -> (i64, i64, i64) {
    (
        (x / cell).floor() as i64,
        (y / cell).floor() as i64,
        (z / cell).floor() as i64,
    )
}
