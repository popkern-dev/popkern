//! Which wedge an imprinted edge enters at a pinch use (PK-1129).
//!
//! A face touching itself at a vertex has its boundary pass through the
//! vertex twice (`join_loops_at_vertex`), so two of its fins leave the
//! vertex, and its interior meets the vertex in two wedges: each from
//! a leaving fin counter-clockwise, about the face's folded normal, to
//! the fin arriving before it (a loop keeps the face on its left). An
//! edge imprinted from the vertex enters one wedge, and must be spliced
//! in before that wedge's leaving fin; the other splice links the
//! wrong sectors and the face's loops no longer bound its interior.
//!
//! The wedge is read off chords, not tangents: a pinch made by a
//! tangency (a dimple's rim touching a plate's edge) has both wedges
//! cusps, every fin and the curve leaving along one line, and the
//! curves are told apart only at second order. Each is cut at one
//! chord distance `s` from the vertex, a hundredth of the shortest
//! chord to a curve's middle, and the cuts' directions in the tangent
//! plane order them. The choice must hold again at `s / 10`, and an
//! edge entering no single wedge, or entering it within a nanoradian
//! of a bounding fin, declines: an osculating curve is not decided.

use super::EulerError;
use crate::geometry::curve::{CurveEval, CurveOnEdge};
use crate::geometry::uv_projection::project_point_to_uv;
use crate::geometry::{Point3, PortableMath, Vec3};
use crate::topology::{BRep, FaceHandle, HalfEdgeHandle, VertexHandle};

/// The chord distance, as a fraction of the shortest chord to a middle.
const CUT_FRACTION: f64 = 1e-2;

/// Angles this close to a bounding fin are on it, not inside.
const ANGLE_MARGIN: f64 = 1e-9;

/// A curve traversed from parameter `from` (at the vertex) towards `to`.
struct Run<'a> {
    curve: &'a CurveOnEdge,
    from: f64,
    to: f64,
}

impl Run<'_> {
    /// The distance from `p` to the run's middle.
    fn reach(&self, p: &Point3) -> f64 {
        self.curve
            .geometry
            .point_at(0.5 * (self.from + self.to))
            .distance_to(p)
    }

    /// The first point of the run at distance `s` from `p`: the run is
    /// scanned outward from the vertex in halvings, then bisected.
    fn cut(&self, p: &Point3, s: f64) -> Option<Point3> {
        let at = |f: f64| {
            self.curve
                .geometry
                .point_at(self.from + f * 0.5 * (self.to - self.from))
        };
        let (mut lo, mut hi) = (0.0, None);
        for k in (0..=48).rev() {
            let f = 0.5_f64.powi(k);
            if at(f).distance_to(p) >= s {
                hi = Some(f);
                break;
            }
            lo = f;
        }
        let mut hi = hi?;
        for _ in 0..80 {
            let mid = 0.5 * (lo + hi);
            if at(mid).distance_to(p) >= s {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        Some(at(hi))
    }
}

impl BRep {
    /// Of `leaving`, the fins of `face`'s loops leaving `v` (more than
    /// one: a pinch use), the one before which an edge along `curve`
    /// leaving `v` (at its start when `at_start`) is spliced: the fin
    /// whose wedge holds the curve's departure.
    pub(super) fn pinch_wedge_fin(
        &self,
        face: FaceHandle,
        v: VertexHandle,
        leaving: &[HalfEdgeHandle],
        curve: &CurveOnEdge,
        at_start: bool,
    ) -> Result<HalfEdgeHandle, EulerError> {
        let decline = |reason: &str| EulerError::NotImprintable {
            face,
            reason: reason.into(),
        };
        let undecided = || decline("the edge leaves a pinch vertex into no single wedge");
        let p = self
            .vertices
            .get(v)
            .ok_or(decline("unknown vertex"))?
            .position;
        let surface = self
            .faces
            .get(face)
            .and_then(|f| f.surface.as_ref())
            .ok_or(decline("the face has no surface"))?;
        let n = project_point_to_uv(&surface.geometry, &p)
            .and_then(|(u, w)| surface.folded_normal(u, w))
            .ok_or(decline("no normal at the pinch vertex"))?;
        let new = if at_start {
            Run {
                curve,
                from: curve.t_start,
                to: curve.t_end,
            }
        } else {
            Run {
                curve,
                from: curve.t_end,
                to: curve.t_start,
            }
        };
        // Each candidate's leaving run and the run arriving before it,
        // read backwards from the vertex.
        let mut wedges = Vec::with_capacity(leaving.len());
        for &h in leaving {
            let prev = self.half_edges[h].prev;
            let out = self.fin_run(h, true).ok_or_else(undecided)?;
            let back = self.fin_run(prev, false).ok_or_else(undecided)?;
            wedges.push((h, out, back));
        }
        let shortest = wedges
            .iter()
            .flat_map(|(_, o, b)| [o.reach(&p), b.reach(&p)])
            .fold(new.reach(&p), f64::min);
        if !(shortest.is_finite() && shortest > 0.0) {
            return Err(undecided());
        }
        let choose = |s: f64| -> Option<HalfEdgeHandle> {
            let dir = |run: &Run| {
                let d = run.cut(&p, s)? - p;
                (d - n * d.dot(&n)).normalized()
            };
            let c = dir(&new)?;
            let mut held = None;
            for (h, out, back) in &wedges {
                let (e1, b) = (dir(out)?, dir(back)?);
                let e2 = n.cross(&e1);
                let angle = |x: &Vec3| {
                    x.dot(&e2)
                        .patan2(x.dot(&e1))
                        .rem_euclid(std::f64::consts::TAU)
                };
                let (ab, ac) = (angle(&b), angle(&c));
                if ac > ANGLE_MARGIN && ac < ab - ANGLE_MARGIN {
                    if held.is_some() {
                        return None;
                    }
                    held = Some(*h);
                }
            }
            held
        };
        let s = CUT_FRACTION * shortest;
        match (choose(s), choose(0.1 * s)) {
            (Some(a), Some(b)) if a == b => Ok(a),
            _ => Err(undecided()),
        }
    }

    /// Fin `h`'s edge curve run from its origin (`from_origin`) or back
    /// from its end, in the fin's sense.
    fn fin_run(&self, h: HalfEdgeHandle, from_origin: bool) -> Option<Run<'_>> {
        let c = self
            .edges
            .get(self.half_edges.get(h)?.edge)?
            .curve
            .as_ref()?;
        let (first, last) = if self.fin_sense(h)? {
            (c.t_start, c.t_end)
        } else {
            (c.t_end, c.t_start)
        };
        Some(if from_origin {
            Run {
                curve: c,
                from: first,
                to: last,
            }
        } else {
            Run {
                curve: c,
                from: last,
                to: first,
            }
        })
    }
}
