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
//! periodic axes unwrapped to the query's branch), voted over points
//! spread round the hole: a hole the split's new edge crosses goes to
//! the side holding most of it. One probe (its first edge's middle)
//! read a fillet's foot circle running 0.01 past a plate edge round's
//! rail on the round's side, half the time by where the circle started,
//! and the foot was bound to the round (PK-1189). A hole that neither
//! or both sides claim as often is left where it is — a contradiction
//! the caller's postconditions must surface, never silently resolved.

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
        // Each side's outer boundary, charted once per surface.
        let mut outlines: Vec<(SurfaceGeometry, Option<Outline>, Option<Outline>)> = Vec::new();
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
            let index = match outlines.iter().position(|(s, _, _)| *s == surface) {
                Some(i) => i,
                None => {
                    let a = self.outer_loop_outline(face_a, &surface);
                    let b = self.outer_loop_outline(face_b, &surface);
                    outlines.push((surface.clone(), a, b));
                    outlines.len() - 1
                }
            };
            let (_, outline_a, outline_b) = &outlines[index];
            let (Some(outline_a), Some(outline_b)) = (outline_a, outline_b) else {
                continue;
            };
            let (mut votes_a, mut votes_b) = (0usize, 0usize);
            for query in self.loop_probe_uvs(hole, &surface) {
                match (outline_a.contains(query), outline_b.contains(query)) {
                    (true, false) => votes_a += 1,
                    (false, true) => votes_b += 1,
                    _ => {}
                }
            }
            let target = match votes_a.cmp(&votes_b) {
                std::cmp::Ordering::Greater => face_a,
                std::cmp::Ordering::Less => face_b,
                std::cmp::Ordering::Equal => continue,
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

    /// UV points spread round the loop: each half-edge's quarter
    /// stations (pcurve, else its edge curve projected).
    fn loop_probe_uvs(
        &self,
        loop_handle: LoopHandle,
        surface: &SurfaceGeometry,
    ) -> Vec<(f64, f64)> {
        const STATIONS: [f64; 3] = [0.25, 0.5, 0.75];
        let mut out = Vec::new();
        let Some(lp) = self.loops.get(loop_handle) else {
            return out;
        };
        let start = lp.half_edge;
        let mut cursor = start;
        for _ in 0..=self.half_edges.len() {
            let Some(he) = self.half_edges.get(cursor) else {
                break;
            };
            for f in STATIONS {
                let uv = if let Some(pc) = &he.pcurve {
                    let p = pc.point_at(pc.t_start + (pc.t_end - pc.t_start) * f);
                    Some((p.x, p.y))
                } else {
                    let p = self.edges.get(he.edge).and_then(|edge| {
                        edge.curve
                            .as_ref()
                            .map(|c| c.geometry.point_at(c.t_start + (c.t_end - c.t_start) * f))
                    });
                    p.and_then(|p| project_point_to_uv(surface, &p))
                };
                out.extend(uv);
            }
            cursor = he.next;
            if cursor == start {
                break;
            }
        }
        out
    }

    /// The face's outer-loop UV polygon for even–odd containment;
    /// `None` when the boundary cannot be sampled.
    ///
    /// On a periodic chart the boundary is chained CONTINUOUSLY (each
    /// sample on the branch nearest its predecessor — the loop lifted
    /// to the universal cover) and tested against every period image
    /// of the query. Unwrapping each vertex to the query's branch
    /// instead folded a full-turn wall's seams behind the query: no
    /// side of a split cylinder housing claimed its fin footprints
    /// (PK-730). A loop that does not close in the cover (it winds
    /// round the period) has no even–odd interior: `None`.
    fn outer_loop_outline(&self, face: FaceHandle, surface: &SurfaceGeometry) -> Option<Outline> {
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
        Some(Outline {
            polygon,
            period_u,
            period_v,
        })
    }
}

/// An outer loop charted for containment: its polygon lifted to the
/// cover and the chart's periods.
struct Outline {
    polygon: Vec<(f64, f64)>,
    period_u: Option<f64>,
    period_v: Option<f64>,
}

impl Outline {
    /// Even–odd containment of `query`.
    fn contains(&self, query: (f64, f64)) -> bool {
        let polygon = &self.polygon;
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
        for &qu in &images(query.0, lo.0, hi.0, self.period_u) {
            for &qv in &images(query.1, lo.1, hi.1, self.period_v) {
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
        inside
    }
}

#[cfg(test)]
#[path = "brep_inner_loop_rehome_tests.rs"]
mod tests;
