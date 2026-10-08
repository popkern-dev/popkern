//! Inner-loop ownership after a face split.
//!
//! [`BRep::mef`] is purely topological: the face it creates starts with
//! no inner loops and the face it splits keeps every hole it had,
//! whichever side of the new edge the hole lies on. A caller that
//! splits a face carrying holes must therefore re-home each hole to
//! the side whose outer boundary contains it — otherwise a hole can
//! end up on a sliver that is later merged away, and its rim edges
//! are then bound to a face they do not lie on (measured on the
//! flagship boss plate under a G2 fillet: the boundary-ring bridge
//! split the plate top, the boss hole stayed with the corner sliver,
//! the sliver merged into the corner cap, and the boss junction circle
//! reported 1.6 off "its" face — PK-389).
//!
//! Containment is decided in the shared surface's UV chart by the
//! even–odd rule against each side's outer-loop polygon (pcurve
//! samples, falling back to projected edge-curve samples, with
//! periodic axes unwrapped to the query's branch). A hole that neither
//! or both sides claim is left where it is — a contradiction the
//! caller's postconditions must surface, never silently resolved.

use crate::geometry::curve::CurveEval;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::uv_projection::project_point_to_uv;
use crate::topology::{BRep, FaceHandle, LoopHandle};

impl BRep {
    /// Re-home the inner loops of `face_a` and `face_b` (the two sides
    /// of a face split that share one surface) to whichever side's
    /// outer boundary contains them. Returns the loops that moved.
    pub fn rehome_inner_loops_after_split(
        &mut self,
        face_a: FaceHandle,
        face_b: FaceHandle,
    ) -> Vec<LoopHandle> {
        self.rehome_inner_loops_except(face_a, face_b, &[])
    }

    /// [`Self::rehome_inner_loops_after_split`], leaving the loops in
    /// `keep` where they are — the hole a split CUT keeps its owner by
    /// construction, and its probe would sit on the new edge both
    /// sides share.
    pub(crate) fn rehome_inner_loops_except(
        &mut self,
        face_a: FaceHandle,
        face_b: FaceHandle,
        keep: &[LoopHandle],
    ) -> Vec<LoopHandle> {
        let mut moved = Vec::new();
        if face_a == face_b {
            return moved;
        }
        let holes: Vec<(FaceHandle, LoopHandle)> = [face_a, face_b]
            .into_iter()
            .flat_map(|f| {
                self.faces
                    .get(f)
                    .map(|face| {
                        face.inner_loops
                            .iter()
                            .map(move |&l| (f, l))
                            .collect::<Vec<_>>()
                    })
                    .unwrap_or_default()
            })
            .collect();
        for (owner, hole) in holes {
            if keep.contains(&hole) {
                continue;
            }
            let Some(surface) = self
                .faces
                .get(owner)
                .and_then(|f| f.surface.as_ref())
                .map(|s| s.geometry.clone())
            else {
                continue;
            };
            let Some(query) = self.loop_probe_uv(hole, &surface) else {
                continue;
            };
            let in_a = self.outer_loop_contains_uv(face_a, &surface, query);
            let in_b = self.outer_loop_contains_uv(face_b, &surface, query);
            let target = match (in_a, in_b) {
                (Some(true), Some(false)) => face_a,
                (Some(false), Some(true)) => face_b,
                _ => continue,
            };
            if target == owner {
                continue;
            }
            if let Some(from) = self.faces.get_mut(owner) {
                from.inner_loops.retain(|&l| l != hole);
            }
            if let Some(to) = self.faces.get_mut(target) {
                to.inner_loops.push(hole);
            }
            if let Some(l) = self.loops.get_mut(hole) {
                l.face = target;
            }
            moved.push(hole);
        }
        moved
    }

    /// A UV point on the loop (its first half-edge's mid-station).
    fn loop_probe_uv(
        &self,
        loop_handle: LoopHandle,
        surface: &SurfaceGeometry,
    ) -> Option<(f64, f64)> {
        let he = self
            .half_edges
            .get(self.loops.get(loop_handle)?.half_edge)?;
        if let Some(pc) = &he.pcurve {
            let p = pc.point_at(0.5 * (pc.t_start + pc.t_end));
            return Some((p.x, p.y));
        }
        let edge = self.edges.get(he.edge)?;
        let p = match edge.curve.as_ref() {
            Some(c) => c.geometry.point_at(0.5 * (c.t_start + c.t_end)),
            None => self.vertices.get(he.origin)?.position,
        };
        project_point_to_uv(surface, &p)
    }

    /// Even–odd containment of `query` in the face's outer-loop UV
    /// polygon; `None` when the boundary cannot be sampled.
    ///
    /// On a periodic chart the boundary is chained CONTINUOUSLY (each
    /// sample on the branch nearest its predecessor — the loop lifted
    /// to the universal cover) and tested against every period image
    /// of the query. Unwrapping each vertex to the query's branch
    /// instead folded a full-turn wall's seams behind the query: no
    /// side of a split cylinder housing claimed its fin footprints
    /// (PK-730). A loop that does not close in the cover (it winds
    /// round the period) has no even–odd interior: `None`.
    fn outer_loop_contains_uv(
        &self,
        face: FaceHandle,
        surface: &SurfaceGeometry,
        query: (f64, f64),
    ) -> Option<bool> {
        const SAMPLES_PER_HALF_EDGE: usize = 8;
        let start = self.loops.get(self.faces.get(face)?.outer_loop)?.half_edge;
        let (period_u, period_v) = surface.periodicity();
        let branch = |x: f64, near: f64, period: Option<f64>| -> f64 {
            match period {
                Some(p) if p > 0.0 => x + ((near - x) / p).round() * p,
                _ => x,
            }
        };
        let mut polygon: Vec<(f64, f64)> = Vec::new();
        let push = |polygon: &mut Vec<(f64, f64)>, (u, v): (f64, f64)| {
            let q = match polygon.last() {
                Some(&(lu, lv)) => (branch(u, lu, period_u), branch(v, lv, period_v)),
                None => (u, v),
            };
            polygon.push(q);
        };
        let mut cursor = start;
        for _ in 0..=self.half_edges.len() {
            let he = self.half_edges.get(cursor)?;
            if let Some(pc) = &he.pcurve {
                for k in 0..SAMPLES_PER_HALF_EDGE {
                    let t = pc.t_start
                        + (pc.t_end - pc.t_start) * (k as f64 / SAMPLES_PER_HALF_EDGE as f64);
                    let p = pc.point_at(t);
                    push(&mut polygon, (p.x, p.y));
                }
            } else {
                let edge = self.edges.get(he.edge)?;
                let points: Vec<crate::geometry::Point3> = match edge.curve.as_ref() {
                    Some(c) => {
                        let forward = edge.half_edge == cursor;
                        (0..SAMPLES_PER_HALF_EDGE)
                            .map(|k| {
                                let f = k as f64 / SAMPLES_PER_HALF_EDGE as f64;
                                let f = if forward { f } else { 1.0 - f };
                                c.geometry.point_at(c.t_start + (c.t_end - c.t_start) * f)
                            })
                            .collect()
                    }
                    None => vec![self.vertices.get(he.origin)?.position],
                };
                for p in points {
                    push(&mut polygon, project_point_to_uv(surface, &p)?);
                }
            }
            cursor = he.next;
            if cursor == start {
                break;
            }
        }
        if polygon.len() < 3 {
            return None;
        }
        // Closure in the cover: the loop's net shift must be zero.
        let (first, last) = (polygon[0], polygon[polygon.len() - 1]);
        let closes = |a: f64, b: f64, period: Option<f64>| match period {
            Some(p) if p > 0.0 => (branch(a, b, Some(p)) - a).abs() < 0.5 * p,
            _ => true,
        };
        if !closes(first.0, last.0, period_u) || !closes(first.1, last.1, period_v) {
            return None;
        }
        // The lifted polygon is bounded, so an even–odd ray along u is
        // sound in the cover; the query counts once per period image
        // overlapping the polygon (parity then reads containment on
        // the quotient chart — the torus included).
        let (lo, hi) = polygon.iter().fold(
            (
                (f64::INFINITY, f64::INFINITY),
                (f64::NEG_INFINITY, f64::NEG_INFINITY),
            ),
            |(lo, hi), p| {
                (
                    (lo.0.min(p.0), lo.1.min(p.1)),
                    (hi.0.max(p.0), hi.1.max(p.1)),
                )
            },
        );
        let images = |q: f64, lo: f64, hi: f64, period: Option<f64>| -> Vec<f64> {
            match period {
                Some(p) if p > 0.0 => {
                    let k0 = ((lo - q) / p).floor() as i64;
                    let k1 = ((hi - q) / p).ceil() as i64;
                    (k0..=k1).map(|k| q + k as f64 * p).collect()
                }
                _ => vec![q],
            }
        };
        let n = polygon.len();
        let mut inside = false;
        for &qu in &images(query.0, lo.0, hi.0, period_u) {
            for &qv in &images(query.1, lo.1, hi.1, period_v) {
                for i in 0..n {
                    let (u1, v1) = polygon[i];
                    let (u2, v2) = polygon[(i + 1) % n];
                    if (v1 > qv) != (v2 > qv) {
                        let x = u1 + (qv - v1) / (v2 - v1) * (u2 - u1);
                        if qu < x {
                            inside = !inside;
                        }
                    }
                }
            }
        }
        Some(inside)
    }
}
