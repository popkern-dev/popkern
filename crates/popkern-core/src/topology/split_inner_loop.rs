//! Inner-loop trim matching for [`BRep::split_face`] — the hole-rim
//! band split (PK-185/PK-186).
//!
//! A fillet band cut on a face's HOLE boundary has both trim
//! endpoints on one inner loop (interior-rail spur tips at reflex
//! hole corners). `split_face`'s outer-loop matcher cannot see them,
//! and an inner-loop `mef` has an extra correctness obligation the
//! outer-loop case never had: the two arcs of the cut loop are not
//! interchangeable. The band arc (fillet edge + spurs) winds WITH
//! the outer loop — a proper standalone face boundary — while the
//! remainder (the enlarged hole ring) keeps winding against it and
//! must STAY the face's inner loop. `mef` transfers the arc walked
//! from `he_entry` to `he_exit` onto the new face, so entry/exit
//! order decides which piece becomes the band; picking the wrong
//! order inverts the region partition (the hole ring would serve as
//! a face and the band as a hole).

use crate::geometry::Point3;
use crate::geometry::uv_projection::project_point_to_uv;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle};

impl BRep {
    /// NEAR-VERTEX rescue for `split_face` endpoints the strict
    /// coincidence pass missed: an intersection junction reaches the
    /// matcher through TWO independently clipped records (each
    /// face's trim computes the same junction to within clip noise),
    /// so the boundary vertex minted by the neighbor's pre-split can
    /// sit a few µm from THIS trim's endpoint (measured: 2.3e-6 on a
    /// bore rim vs a breached-wall chord — just over the strict
    /// coincidence tolerance, and the whole subtract failed). Takes
    /// the NEAREST loop vertex within a 10× band (or the vertex's
    /// own recorded tolerance), never one already claimed by the
    /// other endpoint. Nearest-within-band cannot mis-pair distinct
    /// model vertices unless features are µm-scale, in which case
    /// the strict pass was already ambiguous.
    pub(crate) fn nearest_loop_vertex_within_band(
        &self,
        loop_start_he: HalfEdgeHandle,
        target: &Point3,
        taken: Option<HalfEdgeHandle>,
        tolerance: &Tolerance,
    ) -> Option<HalfEdgeHandle> {
        let max_iters = self.half_edges.len() + 1;
        let mut best: Option<(f64, HalfEdgeHandle)> = None;
        let mut cursor = loop_start_he;
        let mut iters = 0;
        loop {
            let he_data = self.half_edges.get(cursor)?;
            if Some(cursor) != taken {
                let v = &self.vertices[he_data.origin];
                let band = (tolerance.linear() * 10.0).max(v.tolerance.unwrap_or(0.0));
                let d = v.position.distance_to(target);
                if d <= band && best.map(|(bd, _)| d < bd).unwrap_or(true) {
                    best = Some((d, cursor));
                }
            }
            cursor = he_data.next;
            iters += 1;
            if cursor == loop_start_he || iters >= max_iters {
                break;
            }
        }
        best.map(|(_, he)| he)
    }

    /// Match both trim endpoints on ONE inner loop of `face`,
    /// ordered so the arc from entry to exit is the piece that
    /// winds with the outer loop (the band the new face will own).
    /// Returns `None` when no single inner loop carries both
    /// endpoints, or when the winding of either candidate arc
    /// cannot be measured (no surface projection) — the caller then
    /// keeps its historical precondition error.
    pub(crate) fn match_inner_loop_cut(
        &self,
        face: FaceHandle,
        curve_start: &Point3,
        curve_end: &Point3,
        tolerance: &Tolerance,
    ) -> Option<(HalfEdgeHandle, HalfEdgeHandle)> {
        let face_data = self.faces.get(face)?;
        let surface = &face_data.surface.as_ref()?.geometry;
        let loop_entries = |lh: LoopHandle| -> Option<Vec<(HalfEdgeHandle, Point3)>> {
            let start = self.loops.get(lh)?.half_edge;
            let mut out = Vec::new();
            let mut cur = start;
            for _ in 0..self.half_edges.len() + 1 {
                let he = self.half_edges.get(cur)?;
                out.push((cur, self.vertices.get(he.origin)?.position));
                cur = he.next;
                if cur == start {
                    break;
                }
            }
            Some(out)
        };
        let shoelace = |pts: &[(f64, f64)]| -> f64 {
            let mut a = 0.0;
            for i in 0..pts.len() {
                let (x0, y0) = pts[i];
                let (x1, y1) = pts[(i + 1) % pts.len()];
                a += x0 * y1 - x1 * y0;
            }
            0.5 * a
        };
        let uv_of = |pts: &[Point3]| -> Option<Vec<(f64, f64)>> {
            pts.iter()
                .map(|p| project_point_to_uv(surface, p))
                .collect()
        };
        // Reference winding: the outer loop's signed UV area — from its
        // pcurves when it has them (chart-coherent across a periodic
        // seam: a full cylinder wall's outer loop has only its two seam
        // vertices, which project to one u and measured zero area, so a
        // fin root's reflex-end band could never be cut — PK-730), else
        // from its projected vertices.
        let periods = crate::geometry::chart_periods::geometric_chart_periods(surface);
        let outer_area = match self.loop_pcurve_polyline(face_data.outer_loop, periods) {
            Some(poly) => shoelace(&poly),
            None => {
                let outer = loop_entries(face_data.outer_loop)?;
                let outer_pts: Vec<Point3> = outer.iter().map(|(_, p)| *p).collect();
                shoelace(&uv_of(&outer_pts)?)
            }
        };
        if outer_area.abs() < 1e-30 {
            return None;
        }
        for &lh in &face_data.inner_loops {
            let Some(entries) = loop_entries(lh) else {
                continue;
            };
            // Strict pass, then nearest-in-band rescue (the same
            // trust model as the outer matcher: junction vertices
            // minted by a neighbour's pre-split can sit just past
            // strict coincidence).
            let find = |target: &Point3, taken: Option<usize>| -> Option<usize> {
                for (i, (_, p)) in entries.iter().enumerate() {
                    if Some(i) != taken && tolerance.points_coincident(p, target) {
                        return Some(i);
                    }
                }
                let band = tolerance.linear().max(1e-9) * 10.0;
                let mut best: Option<(f64, usize)> = None;
                for (i, (_, p)) in entries.iter().enumerate() {
                    if Some(i) == taken {
                        continue;
                    }
                    let d = p.distance_to(target);
                    if d <= band && best.map(|(bd, _)| d < bd).unwrap_or(true) {
                        best = Some((d, i));
                    }
                }
                best.map(|(_, i)| i)
            };
            let Some(ie) = find(curve_start, None) else {
                continue;
            };
            let Some(ix) = find(curve_end, Some(ie)) else {
                continue;
            };
            // Winding of the entry → exit arc (plus the closing
            // chord): must match the outer loop's sign.
            let n = entries.len();
            let arc_pts = |from: usize, to: usize| -> Vec<Point3> {
                let mut pts = Vec::new();
                let mut i = from;
                loop {
                    pts.push(entries[i].1);
                    if i == to {
                        break;
                    }
                    i = (i + 1) % n;
                }
                pts
            };
            let forward = arc_pts(ie, ix);
            if forward.len() < 3 {
                // A two-vertex arc has no measurable winding; the
                // complementary arc carries the whole ring — it must
                // stay the hole, so the degenerate arc is the band.
                return Some((entries[ie].0, entries[ix].0));
            }
            // Measured in the same chart as the outer loop: pcurve
            // samples of the arc's half-edges when the outer used its
            // pcurves (mixing pcurve and projected charts inverted the
            // band choice on a cylinder wall — PK-730).
            let arc_hes: Vec<HalfEdgeHandle> = {
                let mut v = Vec::new();
                let mut i = ie;
                while i != ix {
                    v.push(entries[i].0);
                    i = (i + 1) % n;
                }
                v
            };
            let area = match (
                self.loop_pcurve_polyline(face_data.outer_loop, periods),
                self.half_edges_pcurve_polyline(&arc_hes, periods),
            ) {
                (Some(_), Some(poly)) => shoelace(&poly),
                _ => shoelace(&uv_of(&forward)?),
            };
            if area * outer_area >= 0.0 {
                return Some((entries[ie].0, entries[ix].0));
            }
            return Some((entries[ix].0, entries[ie].0));
        }
        None
    }

    /// The loop's chart polyline from its half-edges' pcurves (8
    /// samples each), each oriented to continue the previous one;
    /// `None` when any half-edge lacks a pcurve.
    fn loop_pcurve_polyline(
        &self,
        lh: LoopHandle,
        periods: (Option<f64>, Option<f64>),
    ) -> Option<Vec<(f64, f64)>> {
        let start = self.loops.get(lh)?.half_edge;
        let mut hes = Vec::new();
        let mut cur = start;
        for _ in 0..self.half_edges.len() + 1 {
            hes.push(cur);
            cur = self.half_edges.get(cur)?.next;
            if cur == start {
                break;
            }
        }
        self.half_edges_pcurve_polyline(&hes, periods)
    }

    /// The chart polyline through consecutive half-edges' pcurves (8
    /// samples each), each piece oriented and period-shifted to
    /// continue the previous one; `None` when any lacks a pcurve. A
    /// loop crossing a periodic wall's seam carries pcurves authored
    /// in both charts (u near 0 and near 2π); chained unshifted, the
    /// full-period jump corrupts the shoelace winding.
    fn half_edges_pcurve_polyline(
        &self,
        hes: &[HalfEdgeHandle],
        periods: (Option<f64>, Option<f64>),
    ) -> Option<Vec<(f64, f64)>> {
        let mut pieces: Vec<Vec<(f64, f64)>> = Vec::new();
        for &cur in hes {
            let he = self.half_edges.get(cur)?;
            let pc = he.pcurve.as_ref()?;
            let (t0, t1) = pc.parameter_range();
            pieces.push(
                (0..=8)
                    .map(|k| {
                        let p = pc.point_at(t0 + (t1 - t0) * f64::from(k) / 8.0);
                        (p.x, p.y)
                    })
                    .collect(),
            );
        }
        // The period image of `q` nearest `target`.
        let unwrap = |q: (f64, f64), target: (f64, f64)| -> (f64, f64) {
            let near = |x: f64, t: f64, period: Option<f64>| match period {
                Some(p) if p > 0.0 => x + ((t - x) / p).round() * p,
                _ => x,
            };
            (
                near(q.0, target.0, periods.0),
                near(q.1, target.1, periods.1),
            )
        };
        let d = |a: (f64, f64), b: (f64, f64)| ((a.0 - b.0).powi(2) + (a.1 - b.1).powi(2)).sqrt();
        // Chart distance from `q` to `target` modulo the periods.
        let dp = |q: (f64, f64), target: (f64, f64)| d(unwrap(q, target), target);
        // Orient the first piece toward the second, then chain.
        if pieces.len() >= 2 {
            let (f, s) = (&pieces[0], &pieces[1]);
            let tail = f[f.len() - 1];
            let head = f[0];
            let near = |q: (f64, f64)| dp(s[0], q).min(dp(s[s.len() - 1], q));
            if near(head) < near(tail) {
                pieces[0].reverse();
            }
        }
        let mut poly: Vec<(f64, f64)> = Vec::new();
        for mut piece in pieces {
            if let Some(&last) = poly.last() {
                if dp(piece[piece.len() - 1], last) < dp(piece[0], last) {
                    piece.reverse();
                }
                // One shift for the whole piece keeps it continuous.
                let first = piece[0];
                let shifted = unwrap(first, last);
                let (du, dv) = (shifted.0 - first.0, shifted.1 - first.1);
                for q in &mut piece {
                    q.0 += du;
                    q.1 += dv;
                }
            }
            poly.extend(piece);
        }
        Some(poly)
    }
}
