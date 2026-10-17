//! Pole-side witness for faces the winding audit cannot measure
//! (PK-882).
//!
//! A loop that encircles or touches a chart pole (a cone's apex, a
//! sphere's pole) has no meaningful UV winding, so tier 1 of the audit
//! reports it indeterminate, and the tangency witness resolves it only
//! across a G1 boundary. A cone's lateral face meets its cap at a
//! crease: `BRep::cone` wired its base rim backwards for both faces for
//! its whole life with every validator green (PK-881).
//!
//! The interior-left rule needs no winding: walking a face's own
//! half-edge with the flag-folded normal up, the face lies on the left.
//! A face whose loop TOUCHES a pole spans the chart from that pole to
//! its boundary, so at the boundary point farthest (in the pole's chart
//! coordinate) from the pole, the left step must move TOWARD the pole.
//! A loop that encircles a pole without touching it leaves which pole
//! the face holds undecided, and a loop touching two poles (a full
//! sphere's seam spine) has no far side: both stay indeterminate.

use crate::geometry::curve::CurveEval;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::topology::{BRep, FaceHandle, HalfEdgeHandle};

/// Interior samples per boundary edge searched for the far point.
const EDGE_STATIONS: usize = 8;

impl BRep {
    /// `Some(honest)` from the pole-side rule (module docs) for a face
    /// whose outer loop touches exactly one chart pole, else `None`.
    pub(super) fn pole_side_witness(
        &self,
        face: FaceHandle,
        surface: &SurfaceGeometry,
        orientation: bool,
    ) -> Option<bool> {
        let start = self.loops.get(self.faces.get(face)?.outer_loop)?.half_edge;
        let members = self.loop_half_edges(start)?;
        // Pole samples: the loop's uv samples where the chart normal
        // collapses (the audit's own pole test).
        let mut samples = Vec::new();
        for &he in &members {
            samples.extend(self.half_edge_uv_samples(he, surface).ok()?);
        }
        let mut mags: Vec<f64> = samples
            .iter()
            .map(|&(u, v)| surface.normal_at(u, v).length())
            .collect();
        mags.sort_by(f64::total_cmp);
        let median = *mags.get(mags.len() / 2)?;
        let poles: Vec<f64> = samples
            .iter()
            .filter(|&&(u, v)| surface.normal_at(u, v).length() < median * 1.0e-6)
            .map(|&(_, v)| v)
            .collect();
        let pole_v = *poles.first()?;
        if poles.iter().any(|v| (v - pole_v).abs() > 1.0e-6) {
            return None; // two poles: no far side
        }
        // The boundary point farthest from the pole in v, on a
        // non-seam use with a curve.
        let mut far: Option<(f64, HalfEdgeHandle, f64)> = None;
        for &he in &members {
            let h = self.half_edges.get(he)?;
            if self
                .loops
                .get(self.half_edges.get(h.twin)?.loop_)
                .map(|l| l.face)
                == Some(face)
            {
                continue; // seam: both sides are this face
            }
            let Some(curve) = self.edges.get(h.edge).and_then(|e| e.curve.as_ref()) else {
                continue;
            };
            for k in 1..EDGE_STATIONS {
                let t =
                    curve.t_start + (curve.t_end - curve.t_start) * k as f64 / EDGE_STATIONS as f64;
                let (_, v) = surface.closest_point(&curve.geometry.point_at(t))?;
                let gap = (v - pole_v).abs();
                if far.is_none_or(|(g, _, _)| gap > g) {
                    far = Some((gap, he, t));
                }
            }
        }
        let (gap, he, t) = far?;
        if gap <= 1.0e-9 {
            return None;
        }
        let h = self.half_edges.get(he)?;
        let edge = self.edges.get(h.edge)?;
        let curve = edge.curve.as_ref()?;
        let p = curve.geometry.point_at(t);
        // The use's 3-D direction: endpoints on an open edge, the
        // canonical bit on a closed one.
        let (a, b) = (
            curve.geometry.point_at(curve.t_start),
            curve.geometry.point_at(curve.t_end),
        );
        let forward = if a.distance_to(&b) <= self.tolerance.linear() {
            edge.half_edge == he
        } else {
            let o = self.vertices.get(h.origin)?.position;
            o.distance_to(&a) <= o.distance_to(&b)
        };
        let dt = (curve.t_end - curve.t_start) * 1.0e-4 * if forward { 1.0 } else { -1.0 };
        let tangent = (curve.geometry.point_at(t + dt) - p).normalized()?;
        let (u, v) = surface.closest_point(&p)?;
        let normal = surface.normal_at(u, v).normalized()? * if orientation { 1.0 } else { -1.0 };
        let reach = (curve.geometry.point_at(t + dt) - p)
            .length()
            .max(self.tolerance.linear());
        let (_, v_left) = surface.closest_point(&(p + normal.cross(&tangent) * reach))?;
        let toward_pole = (v_left - v) * (pole_v - v);
        (toward_pole.abs() > 0.0).then_some(toward_pole > 0.0)
    }

    /// The half-edges of the loop starting at `start`, in order.
    fn loop_half_edges(&self, start: HalfEdgeHandle) -> Option<Vec<HalfEdgeHandle>> {
        let mut out = vec![start];
        let mut cur = self.half_edges.get(start)?.next;
        while cur != start {
            out.push(cur);
            cur = self.half_edges.get(cur)?.next;
            if out.len() > self.half_edges.len() {
                return None;
            }
        }
        Some(out)
    }
}

#[cfg(test)]
#[path = "winding_pole_witness_tests.rs"]
mod tests;
