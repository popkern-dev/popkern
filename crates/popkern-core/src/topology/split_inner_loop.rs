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
        // Reference winding: the outer loop's signed UV area.
        let outer = loop_entries(face_data.outer_loop)?;
        let outer_pts: Vec<Point3> = outer.iter().map(|(_, p)| *p).collect();
        let outer_area = shoelace(&uv_of(&outer_pts)?);
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
            let area = shoelace(&uv_of(&forward)?);
            if area * outer_area >= 0.0 {
                return Some((entries[ie].0, entries[ix].0));
            }
            return Some((entries[ix].0, entries[ie].0));
        }
        None
    }
}
