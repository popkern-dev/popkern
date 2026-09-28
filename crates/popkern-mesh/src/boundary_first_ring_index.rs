//! Spatial index over a face's boundary-ring segments for the lattice
//! filter and the trim filter of the boundary-first builder (PK-589).
//!
//! Each interior lattice point asks three questions of the rings: is it
//! inside the trim region, is it within the clearance band, and is it
//! inside some segment's diametral disk. Answered by scanning every
//! segment, that was O(points × segments). Both grow as the facet bound
//! tightens, so a coil wall at 3° paid about 8× its 6° cost per halving
//! (0.4 s per build at 6°, repeated per encroachment round). Here each
//! query visits only the segments that can decide it. The per-segment
//! arithmetic is the scan's, verbatim, so every answer is identical:
//! crossing parity is order-independent, and the clearance and disk
//! tests are "any segment passes", which candidate completeness
//! preserves.

use super::BoundaryRing;

/// Target cell count of the metric grid (memory bound for thin charts
/// whose long segments cover many cells).
const GRID_CELLS: f64 = 262_144.0;

/// Relative padding on registered extents, so a segment whose computed
/// distance or disk test passes by a rounding hair is still a
/// candidate in the queried cell.
const PAD: f64 = 1e-9;

/// Inclusive cell-index ranges `((x0, x1), (y0, y1))` an item covers.
type CellRange = ((usize, usize), (usize, usize));

/// Compressed cell → segment lists.
struct Buckets {
    offsets: Vec<u32>,
    items: Vec<u32>,
}

impl Buckets {
    /// `ranges[s]` = inclusive cell-index ranges `(lo, hi)` per axis of
    /// item `s` (1-D when `ny == 1`).
    fn build(nx: usize, ny: usize, ranges: &[CellRange]) -> Self {
        let mut counts = vec![0u32; nx * ny + 1];
        for &((x0, x1), (y0, y1)) in ranges {
            for y in y0..=y1 {
                for x in x0..=x1 {
                    counts[y * nx + x + 1] += 1;
                }
            }
        }
        for k in 1..counts.len() {
            counts[k] += counts[k - 1];
        }
        let mut fill = counts.clone();
        let mut items = vec![0u32; counts[nx * ny] as usize];
        for (s, &((x0, x1), (y0, y1))) in ranges.iter().enumerate() {
            for y in y0..=y1 {
                for x in x0..=x1 {
                    let c = y * nx + x;
                    items[fill[c] as usize] = s as u32;
                    fill[c] += 1;
                }
            }
        }
        Self {
            offsets: counts,
            items,
        }
    }

    fn get(&self, cell: usize) -> &[u32] {
        &self.items[self.offsets[cell] as usize..self.offsets[cell + 1] as usize]
    }
}

/// A 1-D or 2-D cell lattice over `[lo, lo + n·size)`.
#[derive(Clone, Copy)]
struct Axis {
    lo: f64,
    size: f64,
    n: usize,
}

impl Axis {
    fn new(lo: f64, hi: f64, n: usize) -> Self {
        let n = n.max(1);
        Self {
            lo,
            size: ((hi - lo) / n as f64).max(1e-300),
            n,
        }
    }

    /// Cell of `x`, `None` outside the lattice (or NaN).
    fn cell(&self, x: f64) -> Option<usize> {
        let k = ((x - self.lo) / self.size).floor();
        (k >= 0.0 && k < self.n as f64).then_some(k as usize)
    }

    /// Inclusive cell range covering `[a, b]`, clamped to the lattice.
    fn span(&self, a: f64, b: f64) -> (usize, usize) {
        let clamp = |x: f64| {
            ((x - self.lo) / self.size)
                .floor()
                .clamp(0.0, (self.n - 1) as f64) as usize
        };
        (clamp(a), clamp(b))
    }
}

/// Segment-bucketed view of a face's boundary rings.
pub(super) struct RingIndex {
    /// Raw-UV segments, ring by ring, in ring order.
    segments: Vec<((f64, f64), (f64, f64))>,
    su: f64,
    sv: f64,
    clearance: f64,
    slab: Axis,
    slabs: Buckets,
    gx: Axis,
    gy: Axis,
    near: Buckets,
    disks: Buckets,
}

impl RingIndex {
    /// Index `rings` for queries in the metric `(su, sv)` with the
    /// clearance band `clearance` (metric units).
    pub(super) fn new(rings: &[BoundaryRing], su: f64, sv: f64, clearance: f64) -> Self {
        let mut segments = Vec::new();
        for ring in rings {
            let n = ring.uv.len();
            for i in 0..n {
                segments.push((ring.uv[i], ring.uv[(i + 1) % n]));
            }
        }
        let m = segments.len().max(1);
        let (mut v_lo, mut v_hi) = (f64::INFINITY, f64::NEG_INFINITY);
        for &(a, b) in &segments {
            v_lo = v_lo.min(a.1.min(b.1));
            v_hi = v_hi.max(a.1.max(b.1));
        }
        if !(v_lo.is_finite() && v_hi.is_finite()) {
            (v_lo, v_hi) = (0.0, 1.0);
        }
        let slab = Axis::new(
            v_lo,
            v_hi + PAD * (1.0 + v_hi.abs()),
            (m / 2).clamp(1, 4096),
        );
        let slab_ranges: Vec<_> = segments
            .iter()
            .map(|&(a, b)| ((0, 0), slab.span(a.1.min(b.1), a.1.max(b.1))))
            .collect();
        let slabs = Buckets::build(1, slab.n, &slab_ranges);

        // Metric grid over every segment's clearance band and disk.
        let metric = |p: (f64, f64)| (p.0 * su, p.1 * sv);
        let disk = |a: (f64, f64), b: (f64, f64)| {
            let (a, b) = (metric(a), metric(b));
            let c = (0.5 * (a.0 + b.0), 0.5 * (a.1 + b.1));
            let r2 = 0.25 * ((b.0 - a.0).powi(2) + (b.1 - a.1).powi(2));
            (c, (r2 * 1.0201).sqrt())
        };
        let mut boxes = Vec::with_capacity(segments.len());
        let (mut x0, mut x1, mut y0, mut y1) = (
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
        );
        for &(a, b) in &segments {
            let (pa, pb) = (metric(a), metric(b));
            let reach = clearance.max(0.0);
            let near = (
                pa.0.min(pb.0) - reach,
                pa.0.max(pb.0) + reach,
                pa.1.min(pb.1) - reach,
                pa.1.max(pb.1) + reach,
            );
            let (c, r) = disk(a, b);
            let round = (c.0 - r, c.0 + r, c.1 - r, c.1 + r);
            for bx in [near, round] {
                x0 = x0.min(bx.0);
                x1 = x1.max(bx.1);
                y0 = y0.min(bx.2);
                y1 = y1.max(bx.3);
            }
            boxes.push((near, round));
        }
        if !(x0.is_finite() && x1.is_finite() && y0.is_finite() && y1.is_finite()) {
            (x0, x1, y0, y1) = (0.0, 1.0, 0.0, 1.0);
        }
        let pad = PAD * (1.0 + x1.abs().max(x0.abs()).max(y1.abs()).max(y0.abs()));
        let (x0, x1, y0, y1) = (x0 - pad, x1 + pad, y0 - pad, y1 + pad);
        let side = ((x1 - x0) * (y1 - y0) / GRID_CELLS).sqrt().max(1e-300);
        let nx = (((x1 - x0) / side).ceil() as usize).clamp(1, 4096);
        let ny = (((y1 - y0) / side).ceil() as usize).clamp(1, 4096);
        let (gx, gy) = (Axis::new(x0, x1, nx), Axis::new(y0, y1, ny));
        let range = |b: (f64, f64, f64, f64)| {
            let grow = PAD * (1.0 + (b.1 - b.0).abs() + (b.3 - b.2).abs()) + pad;
            (
                gx.span(b.0 - grow, b.1 + grow),
                gy.span(b.2 - grow, b.3 + grow),
            )
        };
        let near_ranges: Vec<_> = boxes.iter().map(|&(n, _)| range(n)).collect();
        let disk_ranges: Vec<_> = boxes.iter().map(|&(_, d)| range(d)).collect();
        Self {
            near: Buckets::build(gx.n, gy.n, &near_ranges),
            disks: Buckets::build(gx.n, gy.n, &disk_ranges),
            segments,
            su,
            sv,
            clearance,
            slab,
            slabs,
            gx,
            gy,
        }
    }

    /// Even-odd containment in the trim region (the ring polygons),
    /// identical to scanning every segment.
    pub(super) fn inside(&self, p: (f64, f64)) -> bool {
        let Some(k) = self.slab.cell(p.1) else {
            return false;
        };
        let mut inside = false;
        for &s in self.slabs.get(k) {
            let ((x1, y1), (x2, y2)) = self.segments[s as usize];
            if (y1 > p.1) != (y2 > p.1) {
                let x_int = x1 + (p.1 - y1) / (y2 - y1) * (x2 - x1);
                if p.0 < x_int {
                    inside = !inside;
                }
            }
        }
        inside
    }

    /// Whether the metric distance from `p` to some ring segment is at
    /// most the clearance.
    pub(super) fn within_clearance(&self, p: (f64, f64)) -> bool {
        let (su, sv) = (self.su, self.sv);
        let q = (p.0 * su, p.1 * sv);
        let (Some(cx), Some(cy)) = (self.gx.cell(q.0), self.gy.cell(q.1)) else {
            return false;
        };
        self.near.get(cy * self.gx.n + cx).iter().any(|&s| {
            let (ra, rb) = self.segments[s as usize];
            let a = (ra.0 * su, ra.1 * sv);
            let b = (rb.0 * su, rb.1 * sv);
            let (dx, dy) = (b.0 - a.0, b.1 - a.1);
            let l2 = dx * dx + dy * dy;
            let t = if l2 < 1e-300 {
                0.0
            } else {
                (((q.0 - a.0) * dx + (q.1 - a.1) * dy) / l2).clamp(0.0, 1.0)
            };
            let (qx, qy) = (a.0 + dx * t, a.1 + dy * t);
            ((q.0 - qx).powi(2) + (q.1 - qy).powi(2)).sqrt() <= self.clearance
        })
    }

    /// Whether `(u, v)` lies inside some segment's (slightly inflated)
    /// diametral disk.
    pub(super) fn in_diametral_disk(&self, u: f64, v: f64) -> bool {
        let (su, sv) = (self.su, self.sv);
        let (x, y) = (u * su, v * sv);
        let (Some(cx), Some(cy)) = (self.gx.cell(x), self.gy.cell(y)) else {
            return false;
        };
        self.disks.get(cy * self.gx.n + cx).iter().any(|&s| {
            let (ra, rb) = self.segments[s as usize];
            let a = (ra.0 * su, ra.1 * sv);
            let b = (rb.0 * su, rb.1 * sv);
            let cx = 0.5 * (a.0 + b.0);
            let cy = 0.5 * (a.1 + b.1);
            let r2 = 0.25 * ((b.0 - a.0).powi(2) + (b.1 - a.1).powi(2));
            (x - cx).powi(2) + (y - cy).powi(2) < r2 * 1.0201
        })
    }
}

#[cfg(test)]
mod tests {
    use super::super::{BoundaryRing, dist_to_rings, point_in_region};
    use super::*;
    use popkern_core::geometry::Point3;
    use popkern_core::topology::EdgeHandle;

    fn ring(uv: Vec<(f64, f64)>) -> BoundaryRing {
        let n = uv.len();
        BoundaryRing {
            uv,
            pos: vec![Point3::new(0.0, 0.0, 0.0); n],
            sources: vec![(EdgeHandle::default(), 0.0, 1.0); n],
        }
    }

    /// Every index answer equals the full scan's, on a thin anisotropic
    /// chart (coil-wall-like metric) with a dense outer ring, a circular
    /// hole, a long diagonal segment whose disk covers much of the
    /// chart, and query points on the lattice rows, on ring vertices
    /// and outside the chart.
    #[test]
    fn index_answers_equal_the_segment_scans() {
        let mut outer: Vec<(f64, f64)> = (0..400).map(|i| (i as f64 * 0.01, 0.0)).collect();
        outer.push((4.0, 0.0));
        outer.push((4.0, 0.5));
        outer.extend((0..100).map(|i| (4.0 - i as f64 * 0.04, 0.5 + 0.001 * i as f64)));
        outer.push((0.0, 0.6));
        let hole: Vec<(f64, f64)> = (0..64)
            .map(|k| {
                let t = -(k as f64) * std::f64::consts::TAU / 64.0;
                (2.0 + 0.1 * t.cos(), 0.3 + 0.1 * t.sin())
            })
            .collect();
        let rings = vec![ring(outer), ring(hole)];
        let (su, sv, clearance) = (0.23, 7.27, 0.004);
        let index = RingIndex::new(&rings, su, sv, clearance);
        let scan_disk = |u: f64, v: f64| {
            let (x, y) = (u * su, v * sv);
            rings.iter().any(|r| {
                let m = r.uv.len();
                (0..m).any(|k| {
                    let a = (r.uv[k].0 * su, r.uv[k].1 * sv);
                    let b = (r.uv[(k + 1) % m].0 * su, r.uv[(k + 1) % m].1 * sv);
                    let (cx, cy) = (0.5 * (a.0 + b.0), 0.5 * (a.1 + b.1));
                    let r2 = 0.25 * ((b.0 - a.0).powi(2) + (b.1 - a.1).powi(2));
                    (x - cx).powi(2) + (y - cy).powi(2) < r2 * 1.0201
                })
            })
        };
        let mut queries: Vec<(f64, f64)> = Vec::new();
        for i in -5..=210 {
            for j in -5..=70 {
                queries.push((i as f64 * 0.02, j as f64 * 0.01));
            }
        }
        queries.extend(rings.iter().flat_map(|r| r.uv.iter().copied()));
        let mut hits = [0usize; 3];
        for &(u, v) in &queries {
            let inside = point_in_region((u, v), &rings);
            let near = dist_to_rings((u, v), &rings, su, sv) <= clearance;
            let disk = scan_disk(u, v);
            assert_eq!(index.inside((u, v)), inside, "inside at ({u}, {v})");
            assert_eq!(
                index.within_clearance((u, v)),
                near,
                "clearance at ({u}, {v})"
            );
            assert_eq!(index.in_diametral_disk(u, v), disk, "disk at ({u}, {v})");
            hits[0] += inside as usize;
            hits[1] += near as usize;
            hits[2] += disk as usize;
        }
        // The fixture exercises every answer both ways.
        assert!(hits.iter().all(|&h| h > 0 && h < queries.len()), "{hits:?}");
    }
}
